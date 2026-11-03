from ecd_epitope_scanner import stability


def test_case_pool_reused_without_changing_score_order(monkeypatch, tmp_path):
    made=[]
    real=stability.ThreadPoolExecutor
    def factory(*args, **kwargs):
        pool=real(*args, **kwargs)
        made.append(pool)
        return pool
    monkeypatch.setattr(stability, 'ThreadPoolExecutor', factory)
    monkeypatch.setattr(stability, '_run_foldx_batch', lambda binary,pdb,chunks,*args: [float(chunks[0][0][1])])
    groups=[[('A',n,'A','W')] for n in range(4)]
    with stability.foldx_execution_scope():
        for _ in range(2):
            assert stability.score_many_with_foldx('foldx',tmp_path/'input.pdb',groups,tmp_path,max_workers=2)==[0.,1.,2.,3.]
        assert len(made)==1
    assert stability._CASE_FOLDX_POOL.get() is None
