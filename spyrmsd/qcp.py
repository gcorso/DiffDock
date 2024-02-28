from typing import Tuple

import numpy as np
from scipy import optimize

from .due import Doi, due

due.cite(
    Doi("10.1107/S0108767305015266"),
    path="spyrmsd.qcp",
    description="QCP method",
)


def M_mtx(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    """
    Compute inner product between coordinate matrices.

    Parameters
    ----------
    A: numpy.ndarray
        Coordinates `A`
    B: numpy.ndarray
        Coordinates `B`

    Returns
    -------
    numpy.ndarray
        Inner product of the coordinate matrices `A` and `B`

    Notes
    -----
    The inner product of the coordinate matrices `A` and `B` corresponds to the matrix
    :math:`\\mathbf{M}`. [1]_

    If :math:`S_{xy}` is defined as

    .. math:: S_{xy} = \\sum_i^N x_{B,i} y_{A,i}

    then :math:`\\mathbf{M}` is the :math:`3\\times 3` matrix given by

    .. math::
       \\begin{pmatrix}
            S_{xx} & S_{xy} & S_{xz} \\\\
            S_{yx} & S_{yy} & S_{yz} \\\\
            S_{zx} & S_{zy} & S_{zz} \\\\
       \\end{pmatrix}

    .. [1] D. L. Theobald, *Rapid calculation of RMSDs using a quaternion-based
       characteristic polynomial*, Acta Crys. A **61**, 478-480 (2005).
    """
    return B.T @ A


def K_mtx(M):
    """
    Compute symmetric key matrix.

    Parameters
    ----------
    M : numpy.ndarray
        Inner product between coordinate matrices

    Returns
    -------
    numpy.ndarray
        Symmetric key matrix

    Notes
    -----
    The symmetric key matrix corresponds to the matrix :math:`\\mathbf{K}`. [2]_

    If :math:`S_{xy}` is defined as

    .. math:: S_{xy} = \\sum_i^N x_{B,i} y_{A,i}

    then :math:`\\mathbf{K}` is the :math:`4\\times 4` symmetric matrix given by

    .. math::
       \\begin{pmatrix}
            S_{xx} + S_{yy} + S_{zz} & S_{yz} - S_{zy} & S_{zx} - S_{xz} & S_{xy} - S_{yx} \\\\
            & S_{xx} - S_{yy} - S_{zz} & S_{xy} + S_{yx} & S_{zx} + S_{xz}\\\\
            & & -S_{xx} + S_{yy} - S_{zz} & S_{yz} - S_{zy} \\\\
            &  & & -S_{xx} - S_{yy} + S_{zz} \\\\
       \\end{pmatrix}

    .. [2] D. L. Theobald, *Rapid calculation of RMSDs using a quaternion-based
       characteristic polynomial*, Acta Crys. A **61**, 478-480 (2005).
    """

    assert M.shape == (3, 3)

    S_xx = M[0, 0]
    S_xy = M[0, 1]
    S_xz = M[0, 2]
    S_yx = M[1, 0]
    S_yy = M[1, 1]
    S_yz = M[1, 2]
    S_zx = M[2, 0]
    S_zy = M[2, 1]
    S_zz = M[2, 2]

    # p = plus, m = minus
    S_xx_yy_zz_ppp = S_xx + S_yy + S_zz
    S_yz_zy_pm = S_yz - S_zy
    S_zx_xz_pm = S_zx - S_xz
    S_xy_yx_pm = S_xy - S_yx
    S_xx_yy_zz_pmm = S_xx - S_yy - S_zz
    S_xy_yx_pp = S_xy + S_yx
    S_zx_xz_pp = S_zx + S_xz
    S_xx_yy_zz_mpm = -S_xx + S_yy - S_zz
    S_yz_zy_pp = S_yz + S_zy
    S_xx_yy_zz_mmp = -S_xx - S_yy + S_zz

    return np.array(
        [
            [S_xx_yy_zz_ppp, S_yz_zy_pm, S_zx_xz_pm, S_xy_yx_pm],
            [S_yz_zy_pm, S_xx_yy_zz_pmm, S_xy_yx_pp, S_zx_xz_pp],
            [S_zx_xz_pm, S_xy_yx_pp, S_xx_yy_zz_mpm, S_yz_zy_pp],
            [S_xy_yx_pm, S_zx_xz_pp, S_yz_zy_pp, S_xx_yy_zz_mmp],
        ]
    )


def coefficients(M: np.ndarray, K: np.ndarray) -> Tuple[float, float, float]:
    """
    Compute quaternion polynomial coefficients.

    Parameters
    ----------
    M : numpy.ndarray
        Inner product between coordinate matrices
    K: numpy.ndarray
        Symmetric key matrix

    Returns
    -------
    Tuple[float, float, float]
        Quaternion polynomial coefficients

    Notes
    _____
    Returns only :math:`\\mathbf{M}`- and :math:`\\mathbf{K}`-dependent coefficients
    are returned. :math:`c_4=1` and :math:`c_3=0` are not returned.

    The :math:`\\mathbf{M}`- and :math:`\\mathbf{K}`-dependent quaternion polynomial
    coefficients are given by

    .. math:: c_2 = -2 \\text{ tr}\\left(\\mathbf{M}^T\\mathbf{M}\\right)

    .. math:: c_1 = -8 \\text{ det}(\\mathbf{M})

    .. math:: c_0 = \\text{ det}(\\mathbf{K})

    """

    c2 = -2 * np.trace(M.T @ M)
    c1 = -8 * np.linalg.det(M)  # TODO: Slow?
    c0 = np.linalg.det(K)  # TODO: Slow?

    return c2, c1, c0


def lambda_max(Ga: float, Gb: float, c2: float, c1: float, c0: float) -> float:
    """
    Find largest root of the quaternion polynomial.

    Parameters
    ----------
    Ga: float
        Inner product of structure A
    Gb:
        Inner product of structure B
    c2:
        Coefficient :math:`c_2` of the quaternion polynomial
    c1:
        Coefficient :math:`c_1` of the quaternion polynomial
    c0:
        Coefficient :math:`c_0` of the quaternion polynomial

    Returns
    -------
    float
        Largest root of the quaternion polynomial (:math:`\\lambda_\\text{max}`)
    """

    def P(x):
        """
        Quaternion polynomial
        """
        return x**4 + c2 * x**2 + c1 * x + c0

    def dP(x):
        """
        Fist derivative of the quaternion polynomial
        """
        return 4 * x**3 + 2 * c2 * x + c1

    x0 = (Ga + Gb) * 0.5

    lmax = optimize.newton(P, x0, fprime=dP)

    return lmax


def _lambda_max_eig(K: np.ndarray) -> float:
    """
    Find largest eigenvalue of :math:`K`.

    Parameters
    ----------
    K: np.ndarray
        Symmetric key matrix

    Returns
    -------
    float
        Largest eigenvalue of :math:`K`, :math:`\\lambda_\\text{max}`
    """
    e, _ = np.linalg.eig(K)

    return max(e)


def qcp_rmsd(A: np.ndarray, B: np.ndarray, atol: float = 1e-9) -> float:
    """
    Compute RMSD using the quaternion polynomial method.

    Parameters
    ----------
    A: numpy.ndarray
        Coordinates of structure A
    B: numpy.ndarray
        Coordinates of structure B
    atol: float
        Absolute tolerance parameter (see notes)

    Returns
    -------
    float
        RMSD between structures `A` and `B`

    Raises
    ------
    AssertionError
        If the shape of structures `A` and `B` is different

    Notes
    -----
    If the structures `A` and `B` can be superimposed exactly (i.e. they differ only
    by center-of-mass translations and rotations), we have

        .. math:: G_a + G_b = 2 \\lambda_\\text{max}

    This means that :math:`s = G_a + G_bb - 2 * \\lambda_\\text{max}` can become
    negative because of numerical errors and therefore :math:`\\sqrt{s}` fails.
    In order to avoid this problem, the final RMSD is set to :math:`0`
    if :math:`|s| < atol`.
    """

    assert A.shape == B.shape

    N = A.shape[0]

    Ga = np.trace(A.T @ A)
    Gb = np.trace(B.T @ B)

    M = M_mtx(A, B)
    K = K_mtx(M)

    c2, c1, c0 = coefficients(M, K)

    try:
        # Fast calculation of the largest eigenvalue of K as root of the characteristic
        # polynomial.
        l_max = lambda_max(Ga, Gb, c2, c1, c0)
    except RuntimeError:  # Newton method fails to converge; see GitHub Issue #35
        # Fallback to (slower) explicit calculation of the largest eigenvalue of K
        l_max = _lambda_max_eig(K)

    s = Ga + Gb - 2 * l_max

    if abs(s) < atol:  # Avoid numerical errors when Ga + Gb = 2 * l_max
        rmsd = 0.0
    else:
        rmsd = np.sqrt(s / N)

    return rmsd
