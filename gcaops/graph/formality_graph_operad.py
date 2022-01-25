r"""
Formality graph operad
"""
from .formality_graph_vector import FormalityGraphVector_dict, FormalityGraphModule_dict
from .formality_graph_basis import FormalityGraphOperadBasis

class FormalityGraphOperation_dict(FormalityGraphVector_dict):
    """
    Element of a FormalityGraphOperad (stored as a dictionary).
    """
    def __init__(self, parent, vector):
        """
        Initialize this graph operation.
        """
        if not isinstance(parent, FormalityGraphOperad_dict):
            raise ValueError("parent must be a FormalityGraphOperad_dict")
        super().__init__(parent, vector)

class FormalityGraphOperad_dict(FormalityGraphModule_dict):
    """
    Operad of formality graphs (with elements stored as dictionaries).
    """
    def __init__(self, base_ring):
        """
        Initialize this graph operad.
        """
        graph_basis = FormalityGraphOperadBasis()
        super().__init__(base_ring, graph_basis)
        self.element_class = FormalityGraphOperation_dict

    def __repr__(self):
        """
        Return a string representation of this graph operad.
        """
        return 'Operad of formality graphs over {}'.format(self._base_ring)

def FormalityGraphOperad(base_ring):
    """
    Return the operad of formality graphs over the given ``base_ring``.
    """
    return FormalityGraphOperad_dict(base_ring)

