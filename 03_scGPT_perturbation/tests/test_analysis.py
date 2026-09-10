"""Small regression tests for numerical direction, matching and counterfactual semantics.

Run: python -m unittest discover -s tests -v
These tests require the analysis dependencies, but no model weights or CUDA.
"""
import unittest
import numpy as np
import pandas as pd
import anndata as ad
from scipy.sparse import csr_matrix
from pipeline.src_metrics import calculate_spc, summarize_table, analyze_gene
from pipeline.perturbation import make_counterfactual
from pipeline.cell_metrics import PreparedGene, compute_knn_label_enrichment, cell_pair_key


class AnalysisTests(unittest.TestCase):
    def test_spc_direction_and_zero(self):
        np.testing.assert_allclose(calculate_spc(np.array([3.,1.,0.,1.]),np.array([1.,3.,0.,0.])),[1.,-1.,0.,2.])

    def test_topn_selection_excludes_low_scores(self):
        table=pd.DataFrame({'gene_name':['A','B','C'],'spc_euclidean_mean':[.2,-.8,.6]})
        summary, selected=summarize_table(table,'TARGET',[2])
        self.assertEqual(selected[2]['top1_gene'],'C')
        self.assertEqual(selected[2]['top2_gene'],'A')
        self.assertAlmostEqual(summary['top2_mean_spc_euclidean_mean'],.4)

    def test_knn_excludes_self_but_retains_other_responder(self):
        r=np.array([[0.,1.],[.1,1.]])
        rp=np.array([[9.9,1.],[10.1,1.]])
        nr=np.array([[10.,1.],[10.2,1.]])
        prep=PreparedGene('T',None,'ok','',2,2,2,2,2,2,r,rp,nr)
        result=compute_knn_label_enrichment(prep,1)
        self.assertEqual(result['knn_nr_fraction_shift_mean_euclidean'],1.)
        self.assertEqual(cell_pair_key('sample_with_underscores_2'),'sample_with_underscores')

    def test_counterfactual_changes_only_x_target_and_keeps_groups(self):
        matrix=csr_matrix([[1.,2.],[0.,3.],[4.,5.]])
        obj=ad.AnnData(matrix,obs=pd.DataFrame({'Response_Final':['Good','Good','Poor']},index=['a_x','b_y','c_z']),var=pd.DataFrame(index=['T','G']))
        obj.layers['data']=matrix.copy()
        combined=make_counterfactual(obj,'T')
        self.assertEqual(combined.obs.group.tolist(),['nonresponder','responder','responder','responder_pert','responder_pert'])
        np.testing.assert_array_equal(combined.X.toarray(),[[4,5],[1,2],[0,3],[0,2],[0,3]])
        np.testing.assert_array_equal(combined.layers['data'].toarray()[-2:],[[1,2],[0,3]])
        self.assertEqual(combined.obs['T_nonzero'].iloc[-2:].tolist(),[True,False])
        np.testing.assert_array_equal(obj.X.toarray(),matrix.toarray())

    def test_src_drops_cls_and_target_gene(self):
        obj=ad.AnnData(csr_matrix([[0.,1.],[1.,1.],[0.,1.]]),
            obs=pd.DataFrame({'group':['nonresponder','responder','responder_pert'],'T_nonzero':[False,False,True]},index=['n_0','r_1','r_2']),var=pd.DataFrame(index=['T','G']))
        # CLS and T deliberately have unrelated values and must not affect G's score.
        emb=np.array([[[99.,0.],[88.,0.],[0.,1.]],[[1.,99.],[0.,88.],[3.,1.]],[[99.,99.],[88.,88.],[1.,1.]]])
        result,count=analyze_gene(obj,emb,'T')
        self.assertEqual(result.gene_name.tolist(),['G'])
        self.assertEqual(count,1)
        self.assertAlmostEqual(result.spc_euclidean_mean.iloc[0],1.)


if __name__=='__main__':
    unittest.main()
