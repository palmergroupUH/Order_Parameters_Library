# Python standard library:
import os

# Local library:

# Third-party libraries:


KB = float(1 << 10)  # Kilobytes
MB = float(1 << 20)  # Megabytes
GB = float(1 << 30)  # Gigabytes


def status_is_ok(filename):

    # file does not exist !

    if (not os.path.isfile(filename)):

        return False

    # file is empty !

    if (os.stat(filename).st_size == 0):

        return False

    # status is ok

    return True


def get_file_size(filename, units):

    file_size_bytes = os.path.getsize(filename)

    if (units == "MB"):

        return file_size_bytes/MB

    elif (units == "GB"):

        return file_size_bytes/GB


def file_size_is_too_big(filename, threshold_size, units):

    file_size_bytes = os.path.getsize(filename)

    file_size = get_file_size(filename, units)

    if (file_size > threshold_size):

        return True
