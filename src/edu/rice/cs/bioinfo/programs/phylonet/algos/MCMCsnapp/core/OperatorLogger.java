package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.core;

/**
 * Created by wendingqiao on 8/18/15.
 */
public class OperatorLogger {

    public String name;
    public int count;
    public int accept;

    OperatorLogger(String name){
        this.name = name;
        count = 0;
        accept = 0;
    }

    public String toString() {
        double rate = (double)(accept) / (double)(count);
        return "Operation:" + name + "; Used:" + count + "; Accepted:" + accept + " ACrate:" + rate;
    }
}