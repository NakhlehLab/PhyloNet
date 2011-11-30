package edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents;

import javax.swing.text.Element;
import java.util.LinkedList;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 10/3/11
 * Time: 7:05 PM
 * To change this template use File | Settings | File Templates.
 */
public class ParameterIdentSet extends ParameterBase
{
    public final Iterable<String> Elements;

    public final String OriginalSource;

    public ParameterIdentSet(int line, int column, String value) {
        super(line, column);

        OriginalSource = value;

        String[] elements = value.substring(1, value.length()-1).split(",");
        LinkedList<String> llElements = new LinkedList<String>();

        for(int i = 0; i<elements.length; i++)
        {
            llElements.add(elements[i].trim());
        }

        Elements = llElements;

    }

     public ParameterIdentSet(int line, int column, Iterable<String> elements) {
        super(line, column);
         OriginalSource = null;
         Elements = elements;
    }

    public <R, T, E extends Exception> R execute(ParameterAlgo<R, T, E> algo, T input) throws E {
        return algo.forIdentSet(this, input);
    }
}
