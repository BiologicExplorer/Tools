"""
cma_network/__init__.py
Re-exports the public API from cma_network/cma_network.py so the package can
be imported flat:  from cma_network import check_cma_network_membership, ...
"""
from cma_network.cma_network import (
    check_cma_network_membership,
    calculate_cma_score,
    list_cma_network,
)

__all__ = [
    "check_cma_network_membership",
    "calculate_cma_score",
    "list_cma_network",
]
