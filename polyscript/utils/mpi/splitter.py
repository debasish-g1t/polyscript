"""
Work-splitting strategies for MPI distribution.

Provides both round-robin (load-balanced) and chunk-based (contiguous)
splitting so runners can choose the best strategy for their workload.
"""

from __future__ import annotations

from typing import List, TypeVar

T = TypeVar("T")


class WorkSplitter:
    """Distribute a list of work items across *n* ranks.

    Static methods — no state, no instantiation required.
    """

    @staticmethod
    def round_robin(items: List[T], n: int) -> List[List[T]]:
        """Split *items* into *n* chunks by round-robin for load balance.

        >>> WorkSplitter.round_robin([1,2,3,4,5], 3)
        [[1, 4], [2, 5], [3]]
        """
        chunks: List[List[T]] = [[] for _ in range(n)]
        for i, item in enumerate(items):
            chunks[i % n].append(item)
        return chunks

    @staticmethod
    def chunks(items: List[T], chunk_size: int) -> List[List[T]]:
        """Split *items* into contiguous chunks of at most *chunk_size*.

        >>> WorkSplitter.chunks([1,2,3,4,5], 2)
        [[1, 2], [3, 4], [5]]
        """
        return [
            items[i : i + chunk_size] for i in range(0, len(items), chunk_size)
        ]
