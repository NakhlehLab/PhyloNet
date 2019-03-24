package edu.rice.cs.bioinfo.programs.phylonet.algos.supernetwork;

/**
 * Created by IntelliJ IDEA.
 * User: zhujiafan
 * Date: 11/21/18
 * Time: 3:08 PM
 * To change this template use File | Settings | File Templates.
 */
public class SNSolver {
    static public SNSummary Solve(SNProblem problem, SNOptions options) {
        SuperNetwork3.eps = options.eps;
        SuperNetwork3.outgroup = options.outgroup;
        SuperNetwork3.reconcileHeights = options.reconcileHeights;
        SuperNetwork3.trustReticulationTime = options.trustReticulationTime;

        System.out.println("Subproblem size: " + problem.subnetworks.size());

        SuperNetwork3 sn = new SuperNetwork3(problem);
        SNSummary summary = new SNSummary();

        if(options.backboneLeaves != null) sn.SetBackboneLeaves(options.backboneLeaves);
        if(options.buildOrder != null) sn.SetBuildOrder(options.buildOrder);

        summary.inferredNetwork = sn.compute();
        summary.netinfos = new SuperNetwork3(problem).subnetworks_;
        summary.taxaNames = sn.getTaxaNames();

        return summary;
    }
}
