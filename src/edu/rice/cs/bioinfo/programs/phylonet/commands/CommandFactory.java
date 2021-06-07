/*
 * Copyright (c) 2012 Rice University.
 *
 * This file is part of PhyloNet.
 *
 * PhyloNet is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * PhyloNet is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PhyloNet.  If not, see <http://www.gnu.org/licenses/>.
 */

package edu.rice.cs.bioinfo.programs.phylonet.commands;

import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.Parameter;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.SyntaxCommand;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.NetworkNonEmpty;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.Networks;
import edu.rice.cs.bioinfo.library.language.richnewick.reading.RichNewickReader;
import edu.rice.cs.bioinfo.library.programming.Proc3;
import org.reflections.Reflections;

import java.lang.annotation.Annotation;
import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;
import java.util.*;


/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 9/6/11
 * Time: 5:43 PM
 * To change this template use File | Settings | File Templates.
 */
public class CommandFactory {


    public static Command make(SyntaxCommand directive, Map<String,NetworkNonEmpty> sourceIdentToNetwork, Proc3<String, Integer, Integer> errorDetected, RichNewickReader<Networks> rnReader, Random rand)
    {
        return make(directive,sourceIdentToNetwork,null,errorDetected,rnReader,rand);
    }

    public static Command make(SyntaxCommand directive, Map<String,NetworkNonEmpty> sourceIdentToNetwork, Map<String,String> sourceIdentToDNA, Proc3<String, Integer, Integer> errorDetected, RichNewickReader<Networks> rnReader, Random rand)
    {
        final String lowerCommandName = directive.getName().toLowerCase();

        ArrayList<Parameter> params = new ArrayList<Parameter>();
        for(Parameter p : directive.getParameters())
        {
            params.add(p);
        }

        //Reflections.log = null; // disable logging in reflections.  Else a lot of output to stdout.
        org.reflections.Reflections reflections = new org.reflections.Reflections("edu.rice.cs.bioinfo.programs.phylonet.commands");


        Set<Class<?>> annotated = reflections.getTypesAnnotatedWith(CommandName.class);

        Class matchingCommand = null;
        for(Class c : annotated)
        {
            for(Annotation a : c.getAnnotations())
            {
                if(((CommandName)a).value().toLowerCase().equals(lowerCommandName))
                {
                    if(matchingCommand == null)
                    {
                        matchingCommand = c;
                    }
                    else
                    {
                        throw new IllegalStateException("More than one class in classpath is marked with the command name " + directive.getName());
                    }
                }
            }
        }

        if(matchingCommand == null)
        {
            throw new IllegalArgumentException(String.format("Unknown command name '%s'.", directive.getName()));
        }

        try
        {
            try
            {
                Constructor<Command> constructor = matchingCommand.getConstructor(SyntaxCommand.class, ArrayList.class, Map.class, Proc3.class, RichNewickReader.class);
                return (Command) constructor.newInstance(directive, params, sourceIdentToNetwork, errorDetected, rnReader);
            }
            catch(NoSuchMethodException e)
            {
                try
                {
                    Constructor<Command> constructor = matchingCommand.getConstructor(SyntaxCommand.class, ArrayList.class, Map.class, Proc3.class, RichNewickReader.class, Random.class);
                    return (Command) constructor.newInstance(directive, params, sourceIdentToNetwork, errorDetected, rnReader, rand);
                }
                catch (NoSuchMethodException e2)
                {
                    Constructor<Command> constructor = matchingCommand.getConstructor(SyntaxCommand.class, ArrayList.class, Map.class,Map.class, Proc3.class, RichNewickReader.class);
                    return (Command) constructor.newInstance(directive, params, sourceIdentToNetwork,sourceIdentToDNA, errorDetected, rnReader);
                }
            }
        }
        catch (NoSuchMethodException e)
        {
            throw new RuntimeException(e);
        }
        catch (InvocationTargetException e)
        {
            throw new RuntimeException(e);
        }
        catch (InstantiationException e)
        {
            throw new RuntimeException(e);
        }
        catch (IllegalAccessException e)
        {
            throw new RuntimeException(e);
        }

    }
}
