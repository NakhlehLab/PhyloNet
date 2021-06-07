package edu.rice.cs.bioinfo.programs.phylonet.algos.SNAPPForNetwork.Tests;

import edu.rice.cs.bioinfo.programs.phylonet.algos.SNAPPForNetwork.ComplexCP;
import edu.rice.cs.bioinfo.programs.phylonet.algos.SNAPPForNetwork.FMatrix;
import edu.rice.cs.bioinfo.programs.phylonet.algos.SNAPPForNetwork.MatrixQ;
import edu.rice.cs.bioinfo.programs.phylonet.algos.SNAPPForNetwork.R;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model.BiAllelicGTR;
import org.apache.commons.math3.complex.Complex;
import org.junit.Test;

import static org.junit.Assert.*;

/**
 * Created by IntelliJ IDEA.
 * User: zhujiafan
 * Date: 11/14/17
 * Time: 2:21 PM
 * To change this template use File | Settings | File Templates.
 */
public class MatrixQTest {
    @Test
    public void expQTtx() throws Exception {
        int mx = 5;

        FMatrix fb = new FMatrix(mx, false);
        R r = new R(mx, new int[]{2});
        fb.set(r, 1.0);
        //System.out.println(fb);

        BiAllelicGTR BAGTRModel = new BiAllelicGTR(new double[] {0.5, 0.5}, new double[] {1.0});
        MatrixQ Q = new MatrixQ(BAGTRModel , mx, 0.04);
        FMatrix ft = new FMatrix(mx, false);
        long start = System.currentTimeMillis();

        //for(int i = 0 ; i < 10000 ; i++)
        ft.setMatrix(Q.expQTtx(0.5, fb.getArr(), fb.getmx()));
        System.out.println("Seconds: " + (System.currentTimeMillis()-start)/1000.0);

        double [] correctArr = {
                0.016855143897457722,
                0.014409887360442538,
                3.6122385538858384E-12,
                1.0383736988197993E-12,
                3.5753809283433515E-12,
                -7.411832375667437E-14,
                4.819812882071966E-12,
                2.3683760353143115E-12,
                -8.819295860352043E-14,
                8.391036273204638E-15,
                -1.7262776141097367E-12,
                -3.3829504685735725E-12,
                -5.1442142066106505E-15,
                2.235168334076419E-15,
                -8.233981390539856E-16,
                -2.450758968606073E-14,
                -2.685023385202187E-13,
                -1.8381449272945195E-14,
                -4.940304774317211E-16,
                -3.780311507212809E-18,
        };

        assertEquals(correctArr.length, ft.getArr().length);

        for(int i = 0 ; i < correctArr.length ; i++) {
            assertEquals(correctArr[i], ft.getArr()[i], 1e-6);
        }

        System.out.println(ft);
    }

    @Test
    public void testComplex() throws Exception {
        long start = System.currentTimeMillis();

        //Complex c = new Complex(1.0, 1.0);
        ComplexCP c = new ComplexCP(1.0, 1.0);

        Complex t1 = new Complex(1.0, 0.0);
        Complex t2 = new Complex(0.0, 1.0);
        for(int i = 0 ; i < 10000000 ; i++) {
            //c = c.divide(1.0);
            c._im /= 1.0;
            c._re /= 1.0;
        }
        System.out.println(c);
        System.out.println((System.currentTimeMillis()-start)/1000.0);
    }

}