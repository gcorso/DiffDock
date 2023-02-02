import math

import torch

"""

    Python implementation of several functions for rotation representation and conversion, 
    using the PyTorch library.

    Rotation can be represented in multiple ways, including axis-angle representation, 
    quaternions, and rotation matrices. These functions allow converting between these 
    representations, and they also implement the Kabsch algorithm to calculate the 
    optimal rotation matrix between two sets of 3D points.

    Quaternion to Matrix
    This function takes in a tensor of shape (..., 4) representing quaternions, 
    where the first element is the real part and the next three elements are the 
    imaginary parts. The function then calculates the corresponding rotation matrix, 
    which is a tensor of shape (..., 3, 3).

    Axis-Angle to Quaternion
    This function takes in a tensor of shape (..., 3) representing rotations 
    in axis-angle form, where the magnitude of the vector is the angle turned 
    anticlockwise in radians around the direction of the vector. The function 
    calculates the corresponding quaternion, which is a tensor of shape (..., 4)
    with the real part first.

    Axis-Angle to Matrix
    This function takes in a tensor of shape (..., 3) representing rotations 
    in axis-angle form, and calculates the corresponding rotation matrix, 
    which is a tensor of shape (..., 3, 3). This function uses the axis-angle 
    to quaternion function to first convert the axis-angle representation to 
    quaternions, and then uses the quaternion to matrix function to obtain 
    the rotation matrix.

    Rigid Transform Kabsch 3D (Torch)
    This function calculates the optimal rigid transformation (rotation and translation) 
    between two sets of 3D points, represented as matrices A and B. The input matrices 
    must have the same number of columns (features) and A must have 3 rows (3D points). 
    The function returns the rotation matrix and translation vector of the optimal 
    rigid transformation. The Kabsch algorithm used here is a method to find the 
    optimal rotation matrix that minimizes the sum of squared distances between 
    the corresponding points in A and B after applying the rotation.

    Simple examples:

        Suppose you have a rotation represented by a vector of shape (3,) with values [0.1, 0.2, 0.3]. 
        You can use the axis_angle_to_quaternion function to convert it to a quaternion of shape (4,).

        Suppose you have a quaternion represented by a tensor of shape (1, 4) with values [0.5, 0.2, 0.3, 0.1]. 
        You can use the quaternion_to_matrix function to convert it to a rotation matrix of shape (1, 3, 3).

        Suppose you have two sets of 3D points, A and B, represented by matrices with shapes (3, 5) and (3, 5) respectively. 
        You can use the rigid_transform_Kabsch_3D_torch function to find the optimal rotation matrix and translation vector to transform A to best align with B.
"""

def quaternion_to_matrix(quaternions):
    """
    From https://pytorch3d.readthedocs.io/en/latest/_modules/pytorch3d/transforms/rotation_conversions.html
    Convert rotations given as quaternions to rotation matrices.

    Args:
        quaternions: quaternions with real part first,
            as tensor of shape (..., 4).

    Returns:
        Rotation matrices as tensor of shape (..., 3, 3).
    """
    r, i, j, k = torch.unbind(quaternions, -1)
    two_s = 2.0 / (quaternions * quaternions).sum(-1)

    o = torch.stack(
        (
            1 - two_s * (j * j + k * k),
            two_s * (i * j - k * r),
            two_s * (i * k + j * r),
            two_s * (i * j + k * r),
            1 - two_s * (i * i + k * k),
            two_s * (j * k - i * r),
            two_s * (i * k - j * r),
            two_s * (j * k + i * r),
            1 - two_s * (i * i + j * j),
        ),
        -1,
    )
    return o.reshape(quaternions.shape[:-1] + (3, 3))


def axis_angle_to_quaternion(axis_angle):
    """
    From https://pytorch3d.readthedocs.io/en/latest/_modules/pytorch3d/transforms/rotation_conversions.html
    Convert rotations given as axis/angle to quaternions.

    Args:
        axis_angle: Rotations given as a vector in axis angle form,
            as a tensor of shape (..., 3), where the magnitude is
            the angle turned anticlockwise in radians around the
            vector's direction.

    Returns:
        quaternions with real part first, as tensor of shape (..., 4).
    """
    angles = torch.norm(axis_angle, p=2, dim=-1, keepdim=True)
    half_angles = 0.5 * angles
    eps = 1e-6
    small_angles = angles.abs() < eps
    sin_half_angles_over_angles = torch.empty_like(angles)
    sin_half_angles_over_angles[~small_angles] = (
            torch.sin(half_angles[~small_angles]) / angles[~small_angles]
    )
    # for x small, sin(x/2) is about x/2 - (x/2)^3/6
    # so sin(x/2)/x is about 1/2 - (x*x)/48
    sin_half_angles_over_angles[small_angles] = (
            0.5 - (angles[small_angles] * angles[small_angles]) / 48
    )
    quaternions = torch.cat(
        [torch.cos(half_angles), axis_angle * sin_half_angles_over_angles], dim=-1
    )
    return quaternions


def axis_angle_to_matrix(axis_angle):
    """
    From https://pytorch3d.readthedocs.io/en/latest/_modules/pytorch3d/transforms/rotation_conversions.html
    Convert rotations given as axis/angle to rotation matrices.

    Args:
        axis_angle: Rotations given as a vector in axis angle form,
            as a tensor of shape (..., 3), where the magnitude is
            the angle turned anticlockwise in radians around the
            vector's direction.

    Returns:
        Rotation matrices as tensor of shape (..., 3, 3).
    """
    return quaternion_to_matrix(axis_angle_to_quaternion(axis_angle))


def rigid_transform_Kabsch_3D_torch(A, B):
    # R = 3x3 rotation matrix, t = 3x1 column vector
    # This already takes residue identity into account.

    assert A.shape[1] == B.shape[1]
    num_rows, num_cols = A.shape
    if num_rows != 3:
        raise Exception(f"matrix A is not 3xN, it is {num_rows}x{num_cols}")
    num_rows, num_cols = B.shape
    if num_rows != 3:
        raise Exception(f"matrix B is not 3xN, it is {num_rows}x{num_cols}")


    # find mean column wise: 3 x 1
    centroid_A = torch.mean(A, axis=1, keepdims=True)
    centroid_B = torch.mean(B, axis=1, keepdims=True)

    # subtract mean
    Am = A - centroid_A
    Bm = B - centroid_B

    H = Am @ Bm.T

    # find rotation
    U, S, Vt = torch.linalg.svd(H)

    R = Vt.T @ U.T
    # special reflection case
    if torch.linalg.det(R) < 0:
        # print("det(R) < R, reflection detected!, correcting for it ...")
        SS = torch.diag(torch.tensor([1.,1.,-1.], device=A.device))
        R = (Vt.T @ SS) @ U.T
    assert math.fabs(torch.linalg.det(R) - 1) < 3e-3  # note I had to change this error bound to be higher

    t = -R @ centroid_A + centroid_B
    return R, t
