package edu.rice.cs.bioinfo.programs.phylonet.algos.summarize;
/*
 *@ClassName: NetNodeInfo
 *@Description: This class is the info of nodes needed by labelling method: nonEquivalent
 *@Author: Zhen Cao
 *@Date:  2019-07-17 20:58
 *@Version: 1.0
 */

public class NetNodeInfo {
    private int _label;
    private double _support;

    public NetNodeInfo(int label){
        _label = label;
        _support = 0;
    }

    /**
     * @Description:
     * @Param: support
     * @return: void
     * @Author: Zhen Cao
     * @Date: 2019-07-17
     */
    public void setSupport(double support){
        _support = support;
    }



    /**
     * @Description:
     * @Param: label
     * @return: void
     * @Author: Zhen Cao
     * @Date: 2019-07-17
     */

    public void setLabel(int label){
        _label = label;
    }

    /**
     * @Description:
     * @Param:
     * @return: int
     * @Author: Zhen Cao
     * @Date: 2019-07-17
     */


    public int getLabel(){
        return _label;
    }

    /**
     * @Description:
     * @Param:
     * @return: double
     * @Author: Zhen Cao
     * @Date: 2019-07-17
     */
    public double getSupport(){
        return _support;
    }
}
