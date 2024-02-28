import copy
import math
import os
import shutil
from functools import partial

import wandb
import torch
torch.multiprocessing.set_sharing_strategy('file_system')

import resource
rlimit = resource.getrlimit(resource.RLIMIT_NOFILE)
resource.setrlimit(resource.RLIMIT_NOFILE, (64000, rlimit[1]))

import yaml
from utils.diffusion_utils import t_to_sigma as t_to_sigma_compl, t_to_sigma_individual
from datasets.loader import construct_loader
from utils.parsing import parse_train_args
from utils.training import train_epoch, test_epoch, loss_function, inference_epoch_fix
from utils.utils import save_yaml_file, get_optimizer_and_scheduler, get_model, ExponentialMovingAverage


def train(args, model, optimizer, scheduler, ema_weights, train_loader, val_loader, t_to_sigma, run_dir, val_dataset2):

    loss_fn = partial(loss_function, tr_weight=args.tr_weight, rot_weight=args.rot_weight,
                      tor_weight=args.tor_weight, no_torsion=args.no_torsion, backbone_weight=args.backbone_loss_weight,
                      sidechain_weight=args.sidechain_loss_weight)

    best_val_loss = math.inf
    best_val_inference_value = math.inf if args.inference_earlystop_goal == 'min' else 0
    best_val_secondary_value = math.inf if args.inference_earlystop_goal == 'min' else 0
    best_epoch = 0
    best_val_inference_epoch = 0

    freeze_params = 0
    scheduler_mode = args.inference_earlystop_goal if args.val_inference_freq is not None else 'min'
    if args.scheduler == 'layer_linear_warmup':
        freeze_params = args.warmup_dur * (args.num_conv_layers + 2) - 1
        print("Freezing some parameters until epoch {}".format(freeze_params))

    print("Starting training...")
    for epoch in range(args.n_epochs):
        if epoch % 5 == 0: print("Run name: ", args.run_name)

        if args.scheduler == 'layer_linear_warmup' and (epoch+1) % args.warmup_dur == 0:
            step = (epoch+1) // args.warmup_dur
            if step < args.num_conv_layers + 2:
                print("New unfreezing step")
                optimizer, scheduler = get_optimizer_and_scheduler(args, model, step=step, scheduler_mode=scheduler_mode)
            elif step == args.num_conv_layers + 2:
                print("Unfreezing all parameters")
                optimizer, scheduler = get_optimizer_and_scheduler(args, model, step=step, scheduler_mode=scheduler_mode)
                ema_weights = ExponentialMovingAverage(model.parameters(), decay=args.ema_rate)
        elif args.scheduler == 'linear_warmup' and epoch == args.warmup_dur:
            print("Moving to plateu scheduler")
            optimizer, scheduler = get_optimizer_and_scheduler(args, model, step=1, scheduler_mode=scheduler_mode,
                                                               optimizer=optimizer)

        logs = {}
        train_losses = train_epoch(model, train_loader, optimizer, device, t_to_sigma, loss_fn, ema_weights if epoch > freeze_params else None)
        print("Epoch {}: Training loss {:.4f}  tr {:.4f}   rot {:.4f}   tor {:.4f}   sc {:.4f}  lr {:.4f}"
              .format(epoch, train_losses['loss'], train_losses['tr_loss'], train_losses['rot_loss'],
                      train_losses['tor_loss'], train_losses['sidechain_loss'], optimizer.param_groups[0]['lr']))

        if epoch > freeze_params:
            ema_weights.store(model.parameters())
            if args.use_ema: ema_weights.copy_to(model.parameters()) # load ema parameters into model for running validation and inference
        val_losses = test_epoch(model, val_loader, device, t_to_sigma, loss_fn, args.test_sigma_intervals)
        print("Epoch {}: Validation loss {:.4f}  tr {:.4f}   rot {:.4f}   tor {:.4f}   sc {:.4f}"
              .format(epoch, val_losses['loss'], val_losses['tr_loss'], val_losses['rot_loss'], val_losses['tor_loss'], val_losses['sidechain_loss']))

        if args.val_inference_freq != None and (epoch + 1) % args.val_inference_freq == 0:
            inf_dataset = [val_loader.dataset.get(i) for i in range(min(args.num_inference_complexes, val_loader.dataset.__len__()))]
            inf_metrics = inference_epoch_fix(model, inf_dataset, device, t_to_sigma, args)
            print("Epoch {}: Val inference rmsds_lt2 {:.3f} rmsds_lt5 {:.3f} min_rmsds_lt2 {:.3f} min_rmsds_lt5 {:.3f}"
                  .format(epoch, inf_metrics['rmsds_lt2'], inf_metrics['rmsds_lt5'], inf_metrics['min_rmsds_lt2'], inf_metrics['min_rmsds_lt5']))
            logs.update({'valinf_' + k: v for k, v in inf_metrics.items()}, step=epoch + 1)

        if args.double_val and args.val_inference_freq != None and (epoch + 1) % args.val_inference_freq == 0:
            inf_dataset = [val_dataset2.get(i) for i in range(min(args.num_inference_complexes, val_dataset2.__len__()))]
            inf_metrics2 = inference_epoch_fix(model, inf_dataset, device, t_to_sigma, args)
            print("Epoch {}: Val inference on second validation rmsds_lt2 {:.3f} rmsds_lt5 {:.3f} min_rmsds_lt2 {:.3f} min_rmsds_lt5 {:.3f}"
                  .format(epoch, inf_metrics2['rmsds_lt2'], inf_metrics2['rmsds_lt5'], inf_metrics2['min_rmsds_lt2'], inf_metrics2['min_rmsds_lt5']))
            logs.update({'valinf2_' + k: v for k, v in inf_metrics2.items()}, step=epoch + 1)
            logs.update({'valinfcomb_' + k: (v + inf_metrics[k])/2 for k, v in inf_metrics2.items()}, step=epoch + 1)

        if args.train_inference_freq != None and (epoch + 1) % args.train_inference_freq == 0:
            inf_dataset = [train_loader.dataset.get(i) for i in range(min(min(args.num_inference_complexes, 300), train_loader.dataset.__len__()))]
            inf_metrics = inference_epoch_fix(model, inf_dataset, device, t_to_sigma, args)
            print("Epoch {}: Train inference rmsds_lt2 {:.3f} rmsds_lt5 {:.3f} min_rmsds_lt2 {:.3f} min_rmsds_lt5 {:.3f}"
                  .format(epoch, inf_metrics['rmsds_lt2'], inf_metrics['rmsds_lt5'], inf_metrics['min_rmsds_lt2'], inf_metrics['min_rmsds_lt5']))
            logs.update({'traininf_' + k: v for k, v in inf_metrics.items()}, step=epoch + 1)

        if epoch > freeze_params:
            if not args.use_ema: ema_weights.copy_to(model.parameters())
            ema_state_dict = copy.deepcopy(model.module.state_dict() if device.type == 'cuda' else model.state_dict())
            ema_weights.restore(model.parameters())

        if args.wandb:
            logs.update({'train_' + k: v for k, v in train_losses.items()})
            logs.update({'val_' + k: v for k, v in val_losses.items()})
            logs['current_lr'] = optimizer.param_groups[0]['lr']
            wandb.log(logs, step=epoch + 1)

        state_dict = model.module.state_dict() if device.type == 'cuda' else model.state_dict()
        if args.inference_earlystop_metric in logs.keys() and \
                (args.inference_earlystop_goal == 'min' and logs[args.inference_earlystop_metric] <= best_val_inference_value or
                 args.inference_earlystop_goal == 'max' and logs[args.inference_earlystop_metric] >= best_val_inference_value):
            best_val_inference_value = logs[args.inference_earlystop_metric]
            best_val_inference_epoch = epoch
            torch.save(state_dict, os.path.join(run_dir, 'best_inference_epoch_model.pt'))
            if epoch > freeze_params:
                torch.save(ema_state_dict, os.path.join(run_dir, 'best_ema_inference_epoch_model.pt'))

        if args.inference_secondary_metric is not None and args.inference_secondary_metric in logs.keys() and \
                (args.inference_earlystop_goal == 'min' and logs[args.inference_secondary_metric] <= best_val_secondary_value or
                 args.inference_earlystop_goal == 'max' and logs[args.inference_secondary_metric] >= best_val_secondary_value):
            best_val_secondary_value = logs[args.inference_secondary_metric]
            if epoch > freeze_params:
                torch.save(ema_state_dict, os.path.join(run_dir, 'best_ema_secondary_epoch_model.pt'))

        if val_losses['loss'] <= best_val_loss:
            best_val_loss = val_losses['loss']
            best_epoch = epoch
            torch.save(state_dict, os.path.join(run_dir, 'best_model.pt'))
            if epoch > freeze_params:
                torch.save(ema_state_dict, os.path.join(run_dir, 'best_ema_model.pt'))

        if args.save_model_freq is not None and (epoch + 1) % args.save_model_freq == 0:
            shutil.copyfile(os.path.join(run_dir, 'best_model.pt'),
                            os.path.join(run_dir, f'epoch{epoch+1}_best_model.pt'))

        if scheduler:
            if epoch < freeze_params or (args.scheduler == 'linear_warmup' and epoch < args.warmup_dur):
                scheduler.step()
            elif args.val_inference_freq is not None:
                scheduler.step(best_val_inference_value)
            else:
                scheduler.step(val_losses['loss'])

        torch.save({
            'epoch': epoch,
            'model': state_dict,
            'optimizer': optimizer.state_dict(),
            'ema_weights': ema_weights.state_dict(),
        }, os.path.join(run_dir, 'last_model.pt'))

    print("Best Validation Loss {} on Epoch {}".format(best_val_loss, best_epoch))
    print("Best inference metric {} on Epoch {}".format(best_val_inference_value, best_val_inference_epoch))


def main_function():
    args = parse_train_args()
    if args.config:
        config_dict = yaml.load(args.config, Loader=yaml.FullLoader)
        arg_dict = args.__dict__
        for key, value in config_dict.items():
            if isinstance(value, list):
                for v in value:
                    arg_dict[key].append(v)
            else:
                arg_dict[key] = value
        args.config = args.config.name
    assert (args.inference_earlystop_goal == 'max' or args.inference_earlystop_goal == 'min')
    if args.val_inference_freq is not None and args.scheduler is not None:
        assert (args.scheduler_patience > args.val_inference_freq) # otherwise we will just stop training after args.scheduler_patience epochs
    if args.cudnn_benchmark:
        torch.backends.cudnn.benchmark = True

    if args.wandb:
        wandb.init(
            entity='',
            settings=wandb.Settings(start_method="fork"),
            project=args.project,
            name=args.run_name,
            config=args
        )

    # construct loader
    t_to_sigma = partial(t_to_sigma_compl, args=args)
    train_loader, val_loader, val_dataset2 = construct_loader(args, t_to_sigma, device)
    
    model = get_model(args, device, t_to_sigma=t_to_sigma)
    optimizer, scheduler = get_optimizer_and_scheduler(args, model, scheduler_mode=args.inference_earlystop_goal if args.val_inference_freq is not None else 'min')
    ema_weights = ExponentialMovingAverage(model.parameters(),decay=args.ema_rate)

    if args.restart_dir:
        try:
            dict = torch.load(f'{args.restart_dir}/{args.restart_ckpt}.pt', map_location=torch.device('cpu'))
            if args.restart_lr is not None: dict['optimizer']['param_groups'][0]['lr'] = args.restart_lr
            optimizer.load_state_dict(dict['optimizer'])
            model.module.load_state_dict(dict['model'], strict=True)
            if hasattr(args, 'ema_rate'):
                ema_weights.load_state_dict(dict['ema_weights'], device=device)
            print("Restarting from epoch", dict['epoch'])
        except Exception as e:
            print("Exception", e)
            dict = torch.load(f'{args.restart_dir}/best_model.pt', map_location=torch.device('cpu'))
            model.module.load_state_dict(dict, strict=True)
            print("Due to exception had to take the best epoch and no optimiser")
    elif args.pretrain_dir:
        dict = torch.load(f'{args.pretrain_dir}/{args.pretrain_ckpt}.pt', map_location=torch.device('cpu'))
        model.module.load_state_dict(dict, strict=True)
        print("Using pretrained model", f'{args.pretrain_dir}/{args.pretrain_ckpt}.pt')

    numel = sum([p.numel() for p in model.parameters()])
    print('Model with', numel, 'parameters')

    if args.wandb:
        wandb.log({'numel': numel})

    # record parameters
    run_dir = os.path.join(args.log_dir, args.run_name)
    yaml_file_name = os.path.join(run_dir, 'model_parameters.yml')
    save_yaml_file(yaml_file_name, args.__dict__)
    args.device = device

    train(args, model, optimizer, scheduler, ema_weights, train_loader, val_loader, t_to_sigma, run_dir, val_dataset2)


if __name__ == '__main__':
    device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
    main_function()
