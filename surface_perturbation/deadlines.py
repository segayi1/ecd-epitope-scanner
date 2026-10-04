"""Check wall-time limits inside triangle queries without changing numerics."""
from __future__ import annotations


class _BudgetTree:
    def __init__(self,tree,budget):self.tree=tree;self.budget=budget;self.calls=0
    def __getattr__(self,name):return getattr(self.tree,name)
    def query(self,*args,**kwargs):
        self.budget.remaining()
        return self.tree.query(*args,**kwargs)
    def query_ball_point(self,*args,**kwargs):
        self.calls+=1
        # Check the clock every point; persist checkpoints less frequently.
        if self.budget.elapsed>=self.budget.limit or self.calls%1024==0:
            self.budget.remaining()
        return self.tree.query_ball_point(*args,**kwargs)


def guard_mesh(mesh,budget):
    if budget is None or getattr(mesh,'_deadline_budget',None) is budget:return mesh
    mesh.tree=_BudgetTree(mesh.tree,budget);mesh.ray_tree=_BudgetTree(mesh.ray_tree,budget)
    mesh._deadline_budget=budget
    return mesh
