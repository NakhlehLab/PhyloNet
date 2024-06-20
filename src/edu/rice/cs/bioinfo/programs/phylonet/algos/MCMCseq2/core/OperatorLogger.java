package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq2.core;

/**
 * Created by wendingqiao on 8/18/15.
 */
public class OperatorLogger {

    public String name;
    public int count;
    public int accept;
    public int rejectOP;

    OperatorLogger(String name){
        this.name = name;
        count = 0;
        accept = 0;
        rejectOP = 0;
    }

    public String toString() {
        double rate = (double)(accept) / (double)(count);
        return "Operation:" + name + "; Used:" + count + "; Accepted:" + accept + " ACrate:" + rate  + " RejectOP:" + rejectOP;
    }
}