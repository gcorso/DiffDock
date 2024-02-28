import os

import numpy as np


def format(fname: str) -> str:
    """
    Extract format extension from file name.

    Parameters
    ----------
    fname : str
        File name

    Returns
    -------
    str
        File extension

    Notes
    -----
    The file extension is returned without the `.` character, i.e. for the file
    `path/filename.ext` the string `ext` is returned.

    If a file is compressed, the `.gz` extension is ignored.
    """
    name, ext = os.path.splitext(fname)

    if ext == ".gz":
        _, ext = os.path.splitext(name)

    return ext[1:]  # Remove "."


def molformat(fname: str) -> str:
    """
    Extract an OpenBabel-friendly format from file name.

    Parameters
    ----------
    fname : str
        File name

    Returns
    -------
    str
        File extension in an OpenBabel-friendly format

    Notes
    -----
    File types in OpenBabel do not always correspond to the file extension. This
    function converts the file extension to an OpenBabel file type.

    The following table shows the different conversions performed by this function:

    ========= =========
    Extension File Type
    --------- ---------
    xyz       XYZ
    ========= =========
    """

    ext = format(fname)

    if ext == "xyz":
        # xyz files in OpenBabel are called XYZ
        ext = "XYZ"

    return ext


def deg_to_rad(angle: float) -> float:
    """
    Convert angle in degrees to angle in radians.

    Parameters
    ----------
    angle : float
        Angle (in degrees)

    Returns
    -------
    float
        Angle (in radians)
    """

    return angle * np.pi / 180.0


def rotate(
    v: np.ndarray, angle: float, axis: np.ndarray, units: str = "rad"
) -> np.ndarray:
    """
    Rotate vector.

    Parameters
    ----------
    v: numpy.array
        3D vector to be rotated
    angle : float
        Angle of rotation (in `units`)
    axis : numpy.array
        3D axis of rotation
    units: {"rad", "deg"}
        Units of `angle` (in radians `rad` or degrees `deg`)

    Returns
    -------
    numpy.array
        Rotated vector

    Raises
    ------
    AssertionError
        If the axis of rotation is not a 3D vector
    ValueError
        If `units` is not `rad` or `deg`
    """

    assert len(axis) == 3

    # Ensure rotation axis is normalised
    axis = axis / np.linalg.norm(axis)

    if units.lower() == "rad":
        pass
    elif units.lower() == "deg":
        angle = deg_to_rad(angle)
    else:
        raise ValueError(
            f"Units {units} for angle is not supported. Use 'deg' or 'rad' instead."
        )

    t1 = np.outer(axis, np.inner(axis, v)).T
    t2 = np.cos(angle) * np.cross(np.cross(axis, v), axis)
    t3 = np.sin(angle) * np.cross(axis, v)

    return t1 + t2 + t3


def center_of_geometry(coordinates: np.ndarray) -> np.ndarray:
    """
    Center of geometry.

    Parameters
    ----------
    coordinates: np.ndarray
        Coordinates

    Returns
    -------
    np.ndarray
        Center of geometry
    """

    assert coordinates.shape[1] == 3

    return np.mean(coordinates, axis=0)


def center(coordinates: np.ndarray) -> np.ndarray:
    """
    Center coordinates.

    Parameters
    ----------
    coordinates: np.ndarray
        Coordinates

    Returns
    -------
    np.ndarray
        Centred coordinates
    """

    return coordinates - center_of_geometry(coordinates)
