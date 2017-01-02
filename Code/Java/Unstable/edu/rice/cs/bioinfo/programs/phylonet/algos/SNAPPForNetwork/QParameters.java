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
    public RateModel _rModel;
    public Integer _M;
    public Double _gTheta;

    public QParameters(RateModel rModel, Integer M) {
        _rModel = rModel;
        _M = M;
    }

    public QParameters(RateModel rModel, Integer M, Double theta) {
        _rModel = rModel;
        _M = M;
        _gTheta = theta;
    }
}
