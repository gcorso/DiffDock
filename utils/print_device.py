import os
import torch
# from utils.utils import get_default_device


def get_default_device():
    if torch.cuda.is_available():
        return torch.device('cuda')
    elif torch.backends.mps.is_available():
        # Not all operations implemented in MPS yet
        use_mps = os.environ.get("PYTORCH_ENABLE_MPS_FALLBACK", "0") == "1"
        if use_mps:
            return torch.device('mps')
        else:
            return torch.device('cpu')
    else:
        return torch.device('cpu')


device = get_default_device()
print(f"DiffDock Device: {device}")