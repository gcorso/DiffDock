import torch
from torch import nn

ACTIVATIONS = {
    'relu': nn.ReLU,
    'silu': nn.SiLU
}


def FCBlock(in_dim, hidden_dim, out_dim, layers, dropout, activation='relu'):
    activation = ACTIVATIONS[activation]
    assert layers >= 2
    sequential = [nn.Linear(in_dim, hidden_dim), activation(), nn.Dropout(dropout)]
    for i in range(layers - 2):
        sequential += [nn.Linear(hidden_dim, hidden_dim), activation(), nn.Dropout(dropout)]
    sequential += [nn.Linear(hidden_dim, out_dim)]
    return nn.Sequential(*sequential)


class GaussianSmearing(torch.nn.Module):
    # used to embed the edge distances
    def __init__(self, start=0.0, stop=5.0, num_gaussians=50):
        super().__init__()
        offset = torch.linspace(start, stop, num_gaussians)
        self.coeff = -0.5 / (offset[1] - offset[0]).item() ** 2
        self.register_buffer('offset', offset)

    def forward(self, dist):
        dist = dist.view(-1, 1) - self.offset.view(1, -1)
        return torch.exp(self.coeff * torch.pow(dist, 2))


class AtomEncoder(torch.nn.Module):
    def __init__(self, emb_dim, feature_dims, sigma_embed_dim, lm_embedding_dim=0):
        """

        Parameters
        ----------
        emb_dim
        feature_dims
            first element of feature_dims tuple is a list with the length of each categorical feature,
            and the second is the number of scalar features
        sigma_embed_dim
        lm_embedding_dim
        """
        #
        super(AtomEncoder, self).__init__()
        self.atom_embedding_list = torch.nn.ModuleList()
        self.num_categorical_features = len(feature_dims[0])
        self.additional_features_dim = feature_dims[1] + sigma_embed_dim + lm_embedding_dim
        for i, dim in enumerate(feature_dims[0]):
            emb = torch.nn.Embedding(dim, emb_dim)
            torch.nn.init.xavier_uniform_(emb.weight.data)
            self.atom_embedding_list.append(emb)

        if self.additional_features_dim > 0:
            self.additional_features_embedder = torch.nn.Linear(self.additional_features_dim + emb_dim, emb_dim)

    def forward(self, x):
        x_embedding = 0
        assert x.shape[1] == self.num_categorical_features + self.additional_features_dim
        for i in range(self.num_categorical_features):
            x_embedding += self.atom_embedding_list[i](x[:, i].long())

        if self.additional_features_dim > 0:
            x_embedding = self.additional_features_embedder(torch.cat([x_embedding, x[:, self.num_categorical_features:]], axis=1))
        return x_embedding


class OldAtomEncoder(torch.nn.Module):

    def __init__(self, emb_dim, feature_dims, sigma_embed_dim, lm_embedding_type=None):
        """

        Parameters
        ----------
        emb_dim
        feature_dims
            first element of feature_dims tuple is a list with the length of each categorical feature,
            and the second is the number of scalar features
        sigma_embed_dim
        lm_embedding_type
        """
        #
        super(OldAtomEncoder, self).__init__()
        self.atom_embedding_list = torch.nn.ModuleList()
        self.num_categorical_features = len(feature_dims[0])
        self.num_scalar_features = feature_dims[1] + sigma_embed_dim
        self.lm_embedding_type = lm_embedding_type
        for i, dim in enumerate(feature_dims[0]):
            emb = torch.nn.Embedding(dim, emb_dim)
            torch.nn.init.xavier_uniform_(emb.weight.data)
            self.atom_embedding_list.append(emb)

        if self.num_scalar_features > 0:
            self.linear = torch.nn.Linear(self.num_scalar_features, emb_dim)
        if self.lm_embedding_type is not None:
            if self.lm_embedding_type == 'esm':
                self.lm_embedding_dim = 1280
            else: raise ValueError('LM Embedding type was not correctly determined. LM embedding type: ', self.lm_embedding_type)
            self.lm_embedding_layer = torch.nn.Linear(self.lm_embedding_dim + emb_dim, emb_dim)

    def forward(self, x):
        x_embedding = 0
        if self.lm_embedding_type is not None:
            assert x.shape[1] == self.num_categorical_features + self.num_scalar_features + self.lm_embedding_dim
        else:
            assert x.shape[1] == self.num_categorical_features + self.num_scalar_features
        for i in range(self.num_categorical_features):
            x_embedding += self.atom_embedding_list[i](x[:, i].long())

        if self.num_scalar_features > 0:
            x_embedding += self.linear(x[:, self.num_categorical_features:self.num_categorical_features + self.num_scalar_features])
        if self.lm_embedding_type is not None:
            x_embedding = self.lm_embedding_layer(torch.cat([x_embedding, x[:, -self.lm_embedding_dim:]], axis=1))
        return x_embedding
