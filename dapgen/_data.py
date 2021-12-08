import os


def get_test_data(f):
    """
    Get test data from the data directory.
    """
    data_dir = os.path.join(os.path.dirname(__file__), "../tests/test-data")
    return os.path.join(data_dir, f)