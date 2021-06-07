package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCtopo.summary;

/**
 * Maintains the summarized record for each operation
 *
 * Created by wendingqiao on 8/18/15.
 */
public class OperatorInfo{

    public String name;
    public int count;
    public int accept;

    public OperatorInfo(String name){
        this.name = name;
        count = 0;
        accept = 0;
    }

    public String toString() {
        double rate = (double)(accept) / (double)(count);
        return "Operation:" + name + "; Used:" + count + "; Accepted:" + accept + " ACrate:" + rate;
    }
}