package edu.rice.cs.bioinfo.programs.phylonet.algos.SNAPPForNetwork;

import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model.RateModel;

/**
 * Created by IntelliJ IDEA.
 * User: zhujiafan
 * Date: 12/12/16
 * Time: 1:58 PM
 * To change this template use File | Settings | File Templates.
 */
public class QParameters {
    public final RateModel _rModel;
    public final Integer _M;
    public final Double _gTheta;
    public final MatrixQ _gMatrix;

    public QParameters(RateModel rModel, Integer M) {
        _rModel = rModel;
        _M = M;
        _gTheta = null;
        _gMatrix = null;
    }

    public QParameters(RateModel rModel, Integer M, Double theta) {
        _rModel = rModel;
        _M = M;
        _gTheta = theta;
        if(theta != null)
            _gMatrix = new MatrixQ(rModel, M, theta);
        else
            _gMatrix = null;
    }

}
