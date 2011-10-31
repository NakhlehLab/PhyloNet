package edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents;

import edu.rice.cs.bioinfo.library.language.pyson._1_0.ast.Identifier;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ast.TaxaMap;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ast.TaxaMapEntry;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 10/12/11
 * Time: 2:43 PM
 * To change this template use File | Settings | File Templates.
 */
public class ParameterTaxaMap extends ParameterBase
{
    public final LinkedList<Map.Entry<String, List<String>>> _mappings = new LinkedList<Map.Entry<String, List<String>>>();

    public ParameterTaxaMap(int line, int column, TaxaMap map) {
        super(line, column);

       for(TaxaMapEntry e : map.Entries)
       {
           String key = e.Key.Content;

           List<String> values = new LinkedList<String>();
           for(Identifier value : e.Values)
           {
                values.add(value.Content);
           }
           _mappings.add(new AbstractMap.SimpleEntry<String, List<String>>(key, values));
       }

    }

    public <R, T, E extends Exception> R execute(ParameterAlgo<R, T, E> algo, T input) throws E {

        return algo.forTaxaMap(this, input);

    }
}
