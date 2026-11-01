import json,sys,unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'analysis'))
from timing_model import contrasts
class IdentifiabilityTests(unittest.TestCase):
    def test_variance_nullspace_and_safe_direction(self):
        rng=np.random.RandomState(814)
        x=np.tile([-2.1,-2.0,2.0,2.1],(32,1))+rng.normal(0,.02,(32,4))
        x-=x.mean(axis=1)[:,None];q=contrasts(x);a=x/np.sum(x*x,axis=1)[:,None]
        n=np.ones_like(x)
        for i in range(4):
            for j in range(4):
                if i!=j:n[:,i]*=x[:,i]-x[:,j]
        n/=np.linalg.norm(n,axis=1)[:,None]
        projected=np.einsum('nri,ni,nsi->nrs',q,n,q)
        np.testing.assert_allclose(projected,0,atol=1e-12)
        v=np.exp(rng.normal(0,.2,(32,4)))
        s=np.einsum('nri,ni,nsi->nrs',q,v,q)
        g=np.einsum('nri,ni->nr',q,v*a)
        gn=np.einsum('nri,ni->nr',q,n*a)
        h=np.stack([-gn[:,1],gn[:,0]],axis=1);h/=np.linalg.norm(h,axis=1)[:,None]
        k=h*(np.sum(h*g,axis=1)/np.einsum('nr,nrs,ns->n',h,s,h))[:,None]
        anew=a-np.einsum('nri,nr->ni',q,k)
        for gamma in (-.4,0,.4):
            vg=v+gamma*n
            delta=np.sum(vg*(anew*anew-a*a),axis=1)
            predicted=-np.sum(h*g,axis=1)**2/np.einsum('nr,nrs,ns->n',h,s,h)
            np.testing.assert_allclose(delta,predicted,atol=1e-12)
            self.assertTrue(np.all(delta<=1e-12))
    def test_equal_nll_can_worsen_slope(self):
        x=np.array([-1.,-.92,.92,1.]);p=np.eye(4)-.25-np.outer(x,x)/np.dot(x,x)
        va=np.ones(4);vb=np.array([.1,1.828,.172,1.9])
        np.testing.assert_allclose(p@np.diag(va)@p,p@np.diag(vb)@p,atol=1e-12)
        X=np.stack([np.ones(4),x],axis=1)
        def slope(v):return np.linalg.solve(X.T@np.diag(1/v)@X,X.T@np.diag(1/v))[1]
        a,b=slope(va),slope(vb)
        self.assertGreater(np.sqrt(np.sum(b*b*va)/np.sum(a*a*va)),1.3)
if __name__=='__main__':unittest.main()
