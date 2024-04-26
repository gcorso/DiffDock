import logging
import os
from typing import Tuple, Union, List

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from e3nn import o3
from e3nn.nn import BatchNorm
from e3nn.o3 import TensorProduct, Linear
from torch_scatter import scatter, scatter_mean

from models.layers import FCBlock


def get_irrep_seq(ns, nv, use_second_order_repr, reduce_pseudoscalars):
    if use_second_order_repr:
        irrep_seq = [
            f'{ns}x0e',
            f'{ns}x0e + {nv}x1o + {nv}x2e',
            f'{ns}x0e + {nv}x1o + {nv}x2e + {nv}x1e + {nv}x2o',
            f'{ns}x0e + {nv}x1o + {nv}x2e + {nv}x1e + {nv}x2o + {nv if reduce_pseudoscalars else ns}x0o'
        ]
    else:
        irrep_seq = [
            f'{ns}x0e',
            f'{ns}x0e + {nv}x1o',
            f'{ns}x0e + {nv}x1o + {nv}x1e',
            f'{ns}x0e + {nv}x1o + {nv}x1e + {nv if reduce_pseudoscalars else ns}x0o'
        ]
    return irrep_seq


def irrep_to_size(irrep):
    irreps = irrep.split(' + ')
    size = 0
    for ir in irreps:
        m, (l, p) = ir.split('x')
        size += int(m) * (2 * int(l) + 1)
    return size


class FasterTensorProduct(torch.nn.Module):
    # Implemented by Bowen Jing
    def __init__(self, in_irreps, sh_irreps, out_irreps, **kwargs):
        super().__init__()
        #for ir in in_irreps:
        #    m, (l, p) = ir
        #    assert l in [0, 1], "Higher order in irreps are not supported"
        #for ir in out_irreps:
        #    m, (l, p) = ir
        #    assert l in [0, 1], "Higher order out irreps are not supported"
        assert o3.Irreps(sh_irreps) == o3.Irreps('1x0e+1x1o'), "sh_irreps don't look like 1st order spherical harmonics"
        self.in_irreps = o3.Irreps(in_irreps)
        self.out_irreps = o3.Irreps(out_irreps)

        in_muls = {'0e': 0, '1o': 0, '1e': 0, '0o': 0}
        out_muls = {'0e': 0, '1o': 0, '1e': 0, '0o': 0}
        for (m, ir) in self.in_irreps: in_muls[str(ir)] = m
        for (m, ir) in self.out_irreps: out_muls[str(ir)] = m

        self.weight_shapes = {
            '0e': (in_muls['0e'] + in_muls['1o'], out_muls['0e']),
            '1o': (in_muls['0e'] + in_muls['1o'] + in_muls['1e'], out_muls['1o']),
            '1e': (in_muls['1o'] + in_muls['1e'] + in_muls['0o'], out_muls['1e']),
            '0o': (in_muls['1e'] + in_muls['0o'], out_muls['0o'])
        }
        self.weight_numel = sum(a * b for (a, b) in self.weight_shapes.values())

    def forward(self, in_, sh, weight):
        in_dict, out_dict = {}, {'0e': [], '1o': [], '1e': [], '0o': []}
        for (m, ir), sl in zip(self.in_irreps, self.in_irreps.slices()):
            in_dict[str(ir)] = in_[..., sl]
            if ir[0] == 1: in_dict[str(ir)] = in_dict[str(ir)].reshape(list(in_dict[str(ir)].shape)[:-1] + [-1, 3])
        sh_0e, sh_1o = sh[..., 0], sh[..., 1:]
        if '0e' in in_dict:
            out_dict['0e'].append(in_dict['0e'] * sh_0e.unsqueeze(-1))
            out_dict['1o'].append(in_dict['0e'].unsqueeze(-1) * sh_1o.unsqueeze(-2))
        if '1o' in in_dict:
            out_dict['0e'].append((in_dict['1o'] * sh_1o.unsqueeze(-2)).sum(-1) / np.sqrt(3))
            out_dict['1o'].append(in_dict['1o'] * sh_0e.unsqueeze(-1).unsqueeze(-1))
            out_dict['1e'].append(torch.linalg.cross(in_dict['1o'], sh_1o.unsqueeze(-2), dim=-1) / np.sqrt(2))
        if '1e' in in_dict:
            out_dict['1o'].append(torch.linalg.cross(in_dict['1e'], sh_1o.unsqueeze(-2), dim=-1) / np.sqrt(2))
            out_dict['1e'].append(in_dict['1e'] * sh_0e.unsqueeze(-1).unsqueeze(-1))
            out_dict['0o'].append((in_dict['1e'] * sh_1o.unsqueeze(-2)).sum(-1) / np.sqrt(3))
        if '0o' in in_dict:
            out_dict['1e'].append(in_dict['0o'].unsqueeze(-1) * sh_1o.unsqueeze(-2))
            out_dict['0o'].append(in_dict['0o'] * sh_0e.unsqueeze(-1))

        weight_dict = {}
        start = 0
        for key in self.weight_shapes:
            in_, out = self.weight_shapes[key]
            weight_dict[key] = weight[..., start:start + in_ * out].reshape(
                list(weight.shape)[:-1] + [in_, out]) / np.sqrt(in_)
            start += in_ * out

        if out_dict['0e']:
            out_dict['0e'] = torch.cat(out_dict['0e'], dim=-1)
            out_dict['0e'] = torch.matmul(out_dict['0e'].unsqueeze(-2), weight_dict['0e']).squeeze(-2)

        if out_dict['1o']:
            out_dict['1o'] = torch.cat(out_dict['1o'], dim=-2)
            out_dict['1o'] = (out_dict['1o'].unsqueeze(-2) * weight_dict['1o'].unsqueeze(-1)).sum(-3)
            out_dict['1o'] = out_dict['1o'].reshape(list(out_dict['1o'].shape)[:-2] + [-1])

        if out_dict['1e']:
            out_dict['1e'] = torch.cat(out_dict['1e'], dim=-2)
            out_dict['1e'] = (out_dict['1e'].unsqueeze(-2) * weight_dict['1e'].unsqueeze(-1)).sum(-3)
            out_dict['1e'] = out_dict['1e'].reshape(list(out_dict['1e'].shape)[:-2] + [-1])

        if out_dict['0o']:
            out_dict['0o'] = torch.cat(out_dict['0o'], dim=-1)
            # out_dict['0o'] = (out_dict['0o'].unsqueeze(-1) * weight_dict['0o']).sum(-2)
            out_dict['0o'] = torch.matmul(out_dict['0o'].unsqueeze(-2), weight_dict['0o']).squeeze(-2)

        out = []
        for _, ir in self.out_irreps:
            out.append(out_dict[str(ir)])
        return torch.cat(out, dim=-1)


def tp_scatter_simple(tp, fc_layer, node_attr, edge_index, edge_attr, edge_sh,
                      out_nodes=None, reduce='mean', edge_weight=1.0):
    """
    Perform TensorProduct + scatter operation, aka graph convolution.

    This function is only for edge_groups == 1. For multiple edge groups, and for larger graphs,
    use tp_scatter_multigroup instead.
    """

    assert isinstance(edge_attr, torch.Tensor), \
        "This function is only for a single edge group, so edge_attr must be a tensor and not a list."

    _device = node_attr.device
    _dtype = node_attr.dtype
    edge_src, edge_dst = edge_index
    out_irreps = fc_layer(edge_attr).to(_device).to(_dtype)
    out_irreps.mul_(edge_weight)
    tp = tp(node_attr[edge_dst], edge_sh, out_irreps)
    out_nodes = out_nodes or node_attr.shape[0]
    out = scatter(tp, edge_src, dim=0, dim_size=out_nodes, reduce=reduce)
    return out


def tp_scatter_multigroup(tp: o3.TensorProduct, fc_layer: Union[nn.Module, nn.ModuleList],
                          node_attr: torch.Tensor, edge_index: torch.Tensor,
                          edge_attr_groups: List[torch.Tensor], edge_sh: torch.Tensor,
                          out_nodes=None, reduce='mean', edge_weight=1.0):
    """
    Perform TensorProduct + scatter operation, aka graph convolution.

    To keep the peak memory usage reasonably low, this function does not concatenate the edge_attr_groups.
    Rather, we sum the output of the tensor product for each edge group, and then divide by the number of edges

    Parameters
    ----------
    tp: o3.TensorProduct
    fc_layer: nn.Module, or nn.ModuleList
        If a list, must be the same length as edge_attr_groups
    node_attr: torch.Tensor
    edge_index: torch.Tensor of shape (2, num_edges)
        Indicates the source and destination nodes of each edge
    edge_attr_groups: List[torch.Tensor]
        List of tensors, with shape (X_i, num_edge_attributes). Each tensor is a different group of edge attributes
        X may be different for each tensor, although sum(X_i) must be equal to edge_index.shape[1]
    edge_sh: torch.Tensor
        Spherical harmonics for the edges (see o3.spherical_harmonics)
    out_nodes:
        Number of output nodes
    reduce: str
        'mean' or 'sum'. Reduce function for scatter.
    edge_weight : float or torch.Tensor
        Edge weights. If a tensor, must be the same shape as `edge_index`

    Returns
    -------
    torch.Tensor
        Result of the graph convolution
    """

    assert isinstance(edge_attr_groups, list), "This function is only for a list of edge groups"
    assert reduce in {"mean", "sum"}, "Only 'mean' and 'sum' are supported for reduce"
    # It would be possible to support mul/min/max but that would require more work and more code,
    # so only going to do it if it's needed.

    _device = node_attr.device
    _dtype = node_attr.dtype
    edge_src, edge_dst = edge_index
    edge_attr_lengths = [_edge_attr.shape[0] for _edge_attr in edge_attr_groups]
    total_rows = sum(edge_attr_lengths)
    assert total_rows == edge_index.shape[1], "Sum of edge_attr_groups must be equal to edge_index.shape[1]"
    num_edge_groups = len(edge_attr_groups)
    edge_weight_is_indexable = hasattr(edge_weight, '__getitem__')

    out_nodes = out_nodes or node_attr.shape[0]
    total_output_dim = sum([x.dim for x in tp.irreps_out])
    final_out = torch.zeros((out_nodes, total_output_dim), device=_device, dtype=_dtype)
    div_factors = torch.zeros(out_nodes, device=_device, dtype=_dtype)

    cur_start = 0
    for ii in range(num_edge_groups):
        cur_length = edge_attr_lengths[ii]
        cur_end = cur_start + cur_length
        cur_edge_range = slice(cur_start, cur_end)
        cur_edge_src, cur_edge_dst = edge_src[cur_edge_range], edge_dst[cur_edge_range]

        cur_fc = fc_layer[ii] if isinstance(fc_layer, nn.ModuleList) else fc_layer
        cur_out_irreps = cur_fc(edge_attr_groups[ii])
        if edge_weight_is_indexable:
            cur_out_irreps.mul_(edge_weight[cur_edge_range])
        else:
            cur_out_irreps.mul_(edge_weight)

        summand = tp(node_attr[cur_edge_dst, :], edge_sh[cur_edge_range, :], cur_out_irreps)
        # We take a simple sum, and then add up the count of edges which contribute,
        # so that we can take the mean later.
        final_out += scatter(summand, cur_edge_src, dim=0, dim_size=out_nodes, reduce="sum")
        div_factors += torch.bincount(cur_edge_src, minlength=out_nodes)

        cur_start = cur_end

        del cur_out_irreps, summand

    if reduce == 'mean':
        div_factors = torch.clamp(div_factors, torch.finfo(_dtype).eps)
        final_out = final_out / div_factors[:, None]

    return final_out


class TensorProductConvLayer(torch.nn.Module):
    def __init__(self, in_irreps, sh_irreps, out_irreps, n_edge_features, residual=True, batch_norm=True, dropout=0.0,
                 hidden_features=None, faster=False, edge_groups=1, tp_weights_layers=2, activation='relu', depthwise=False):
        super(TensorProductConvLayer, self).__init__()
        self.in_irreps = in_irreps
        self.out_irreps = out_irreps
        self.sh_irreps = sh_irreps
        self.residual = residual
        self.edge_groups = edge_groups
        self.out_size = irrep_to_size(out_irreps)
        self.depthwise = depthwise
        if hidden_features is None:
            hidden_features = n_edge_features

        if depthwise:
            in_irreps = o3.Irreps(in_irreps)
            sh_irreps = o3.Irreps(sh_irreps)
            out_irreps = o3.Irreps(out_irreps)

            irreps_mid = []
            instructions = []
            for i, (mul, ir_in) in enumerate(in_irreps):
                for j, (_, ir_edge) in enumerate(sh_irreps):
                    for ir_out in ir_in * ir_edge:
                        if ir_out in out_irreps:
                            k = len(irreps_mid)
                            irreps_mid.append((mul, ir_out))
                            instructions.append((i, j, k, "uvu", True))

            # We sort the output irreps of the tensor product so that we can simplify them
            # when they are provided to the second o3.Linear
            irreps_mid = o3.Irreps(irreps_mid)
            irreps_mid, p, _ = irreps_mid.sort()

            # Permute the output indexes of the instructions to match the sorted irreps:
            instructions = [
                (i_in1, i_in2, p[i_out], mode, train)
                for i_in1, i_in2, i_out, mode, train in instructions
            ]

            self.tp = TensorProduct(
                in_irreps,
                sh_irreps,
                irreps_mid,
                instructions,
                shared_weights=False,
                internal_weights=False,
            )

            self.linear_2 = Linear(
                # irreps_mid has uncoallesed irreps because of the uvu instructions,
                # but there's no reason to treat them seperately for the Linear
                # Note that normalization of o3.Linear changes if irreps are coallesed
                # (likely for the better)
                irreps_in=irreps_mid.simplify(),
                irreps_out=out_irreps,
                internal_weights=True,
                shared_weights=True,
            )

        else:
            if faster:
                print("Faster Tensor Product")
                self.tp = FasterTensorProduct(in_irreps, sh_irreps, out_irreps)
            else:
                self.tp = o3.FullyConnectedTensorProduct(in_irreps, sh_irreps, out_irreps, shared_weights=False)

        if edge_groups == 1:
            self.fc = FCBlock(n_edge_features, hidden_features, self.tp.weight_numel, tp_weights_layers, dropout, activation)
        else:
            self.fc = [FCBlock(n_edge_features, hidden_features, self.tp.weight_numel, tp_weights_layers, dropout, activation) for _ in range(edge_groups)]
            self.fc = nn.ModuleList(self.fc)

        self.batch_norm = BatchNorm(out_irreps) if batch_norm else None

    def forward(self, node_attr, edge_index, edge_attr, edge_sh, out_nodes=None, reduce='mean', edge_weight=1.0):
        if edge_index.shape[1] == 0 and node_attr.shape[0] == 0:
            raise ValueError("No edges and no nodes")

        _dtype = node_attr.dtype
        if edge_index.shape[1] == 0:
            out = torch.zeros((node_attr.shape[0], self.out_size), dtype=_dtype, device=node_attr.device)
        else:
            if self.edge_groups == 1:
                out = tp_scatter_simple(self.tp, self.fc, node_attr, edge_index, edge_attr, edge_sh,
                                        out_nodes, reduce, edge_weight)
            else:
                out = tp_scatter_multigroup(self.tp, self.fc, node_attr, edge_index, edge_attr, edge_sh,
                                            out_nodes, reduce, edge_weight)

            if self.depthwise:
                out = self.linear_2(out)

            if self.batch_norm:
                out = self.batch_norm(out)

        if self.residual:
            padded = F.pad(node_attr, (0, out.shape[-1] - node_attr.shape[-1]))
            out = out + padded

        out = out.to(_dtype)
        return out


class OldTensorProductConvLayer(torch.nn.Module):
    def __init__(self, in_irreps, sh_irreps, out_irreps, n_edge_features, residual=True, batch_norm=True, dropout=0.0,
                 hidden_features=None):
        super(OldTensorProductConvLayer, self).__init__()
        self.in_irreps = in_irreps
        self.out_irreps = out_irreps
        self.sh_irreps = sh_irreps
        self.residual = residual
        if hidden_features is None:
            hidden_features = n_edge_features

        self.tp = tp = o3.FullyConnectedTensorProduct(in_irreps, sh_irreps, out_irreps, shared_weights=False)

        self.fc = nn.Sequential(
            nn.Linear(n_edge_features, hidden_features),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_features, tp.weight_numel)
        )
        self.batch_norm = BatchNorm(out_irreps) if batch_norm else None

    def forward(self, node_attr, edge_index, edge_attr, edge_sh, out_nodes=None, reduce='mean', edge_weight=1.0):

        # Break up the edge_attr into chunks to limit the maximum memory usage
        edge_chunk_size = 100_000
        num_edges = edge_attr.shape[0]
        num_chunks = (num_edges // edge_chunk_size) if num_edges % edge_chunk_size == 0 \
            else (num_edges // edge_chunk_size) + 1
        edge_ranges = np.array_split(np.arange(num_edges), num_chunks)
        edge_attr_groups = [edge_attr[cur_range] for cur_range in edge_ranges]

        out = tp_scatter_multigroup(self.tp, self.fc, node_attr, edge_index, edge_attr_groups, edge_sh,
                                    out_nodes, reduce, edge_weight)

        if self.residual:
            padded = F.pad(node_attr, (0, out.shape[-1] - node_attr.shape[-1]))
            out = out + padded

        if self.batch_norm:
            out = self.batch_norm(out)

        out = out.to(node_attr.dtype)
        return out
