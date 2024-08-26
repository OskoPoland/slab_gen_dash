import json
import numpy as np

class NpEncoder(json.JSONEncoder):
    """
    TODO: Add documentation string describing usage of class.
    """

    def default(self, obj):
        """
        TODO: Add documentation string for class method.
        """
        if isinstance(obj, np.integer):
            return int(obj)

        if isinstance(obj, np.floating):
            return float(obj)

        if isinstance(obj, np.ndarray):
            return obj.tolist()

        return json.JSONEncoder.default(self, obj)
