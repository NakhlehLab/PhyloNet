package edu.rice.cs.bioinfo.programs.phylonet.algos.summarize;
/*
 *@ClassName: Utils
 *@Description
 *@Author: Zhen Cao
 *@Date:  2019-07-17 17:35
 *@Version: 1.0
 */


import java.util.*;

public class Utils {




    /*  global variables
    *
    * */
    //internal node: <label, counter>
    public static Map<Integer, Double> _labelCounter = new HashMap<>();

    //internal node: <child set, label>
    public static Map<HashSet<Integer>, Integer> _setLabel = new HashMap<>();
    public static int _counter = 0;

    //Leaf node: <name, label>
    public static Map<String, Integer> _leafLabel = new HashMap<>();



}
