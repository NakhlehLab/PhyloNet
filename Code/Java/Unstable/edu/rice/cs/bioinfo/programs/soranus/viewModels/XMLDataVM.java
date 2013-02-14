package edu.rice.cs.bioinfo.programs.soranus.viewModels;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 2/8/13
 * Time: 2:23 PM
 * To change this template use File | Settings | File Templates.
 */
public class XMLDataVM implements DocumentVM
{
    public final String Content;


    public XMLDataVM(String content) {
        Content = content;
    }

    public <R, E extends Exception> R execute(DocumentVMAlgo<R, E> algo) throws E {
        return algo.forXMLDataVM(this);
    }
}
