/**
 * This file is part of PhyloNet-HMM.
 *
 * Copyright © 2013-2014 Kevin Liu, Jingxuan Dai, Kathy Truong, 
 * Ying Song, Michael H. Kohn, and Luay Nakhleh. <http://bioinfo.cs.rice.edu/>
 * 
 * PhyloNet-HMM is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * PhyloNet-HMM is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package edu.rice.cs.bioinfo.programs.phmm.src.reader;

import java.util.ArrayList;
import java.util.HashMap;

/**
 * This class use the base calculation to translate between letter combinations and integers.
 *
 * @author Jingxuan Dai
 *
 */

public class IntTranslator {

    /**
     * This is a help function that maps each element in the given alphabet to a unique integer starting from 0.
     *
     * @param alphabet: an ArrayList of string, in which each string is an element of the alphabet (sample space).
     * @return A HashMap, in which keys are strings and values are the corresponding integers.
     */

    public static HashMap<String, Integer>  buildMap (ArrayList<String> alphabet) {
        HashMap<String, Integer> aMap = new HashMap<String, Integer> ();
        for (int i = 0; i < alphabet.size(); i++) {
            aMap.put(alphabet.get(i), i);
            // System.out.println(alphabet.get(i) + ", " + i);
        }
        return aMap;
    }

    /**
     * Given a predefined mapping of alphabet elements and integers, this function translates a combination letter sequence to a
     * number.
     *
     * @param combi: The letter combination, which is a string
     * @param aMap: The given mapping of letters and integers, which is a HashMap
     * @return The letter combination's corresponding integer
     */

    public static int letterToInt(String combi, HashMap<String, Integer> aMap) {

        int result = 0;

        for (int i = 0; i < combi.length(); i++) {
            if (aMap.containsKey(combi.substring(i, i+1))) {
                // System.out.println(Math.pow(aMap.size(), combi.length() - 1 - i));
                result += aMap.get(combi.substring(i, i+1)) * (Math.pow(aMap.size(), combi.length() - 1 - i));
            }
            else
                throw new RuntimeException("The letter " + combi.substring(i, i+1) + " in combination " + combi + " does not exist in our alphabet.");
        }

        return result;
    }

    /**
     * Given a predefined mapping, this function translates a number to its corresponding letter sequence.
     *
     * @param numSeq: number of sequences, which is also the number of letters in a letter combination
     * @param num: the number that we want to translate
     * @param aMap : the predefined mapping
     * @return A string, which is the translated letter sequence
     */

    public static String intToLetter (int numSeq, int num, HashMap<String, Integer> aMap) {

        String combi = "";
        String tempNum = "";
        int base = aMap.size();

        while (num >= base) {
            // tempNum.concat(String.valueOf(num%base));
            tempNum = String.valueOf(num%base).concat(tempNum);
            num = num/base;
            // System.out.println(num);
        }

        tempNum = String.valueOf(num).concat(tempNum);

        while (tempNum.length() < numSeq ) {
            String temp = "0";
            tempNum = temp.concat(tempNum);
        }

        // System.out.println(tempNum);

        for (int i = 0; i < tempNum.length(); i++) {
            int digit = Integer.parseInt(tempNum.substring(i, i+1));
            // System.out.println(digit);
            for (String str: aMap.keySet()) {
                if (aMap.get(str) == digit) {
                    // System.out.println(str);
                    combi = combi.concat(str);
                    break;
                }
            }
        }

        return combi;
    }



    /*
    public static void main(String[] args) {
        // TODO Auto-generated method stub

        ArrayList<String> alphabet = new ArrayList<String> ();

        alphabet.add("A");
        alphabet.add("T");
        alphabet.add("G");
        alphabet.add("C");

        HashMap<String, Integer> aMap = buildMap(alphabet);
        System.out.println(intToLetter(4, 0, aMap));
        // String combi = "AAAAA";
        // System.out.println(letterToInt(combi, aMap));

    }
    */

}
