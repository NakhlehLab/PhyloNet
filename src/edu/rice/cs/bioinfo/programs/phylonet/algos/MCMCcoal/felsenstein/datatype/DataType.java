package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.felsenstein.datatype;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Created by wendingqiao on 5/3/16.
 */
public interface DataType {
    final static public char GAP_CHAR = '-';
    final static public char MISSING_CHAR = '?';

    /**
     * gets number of states.
     */
    int getStateCount();

    /**
     * Convert a sequence represented by a string into a sequence of data type integers
     */
    List<Integer> stringToState(String sequence);

    /**
     * Convert a list of states into a sequence represented by a string.
     */
    String stateToString(List<Integer> states);

    /**
     * Convert an array of states into a sequence represented by a string.
     */
    String stateToString(int[] states);

    /**
     * returns an array of length getStateCount().
     */
    boolean[] getStateSet(int state);

    /**
     * returns an array with all non-ambiguous states.
     */
    int[] getStatesForCode(int state);

    /**
     * returns whether the state is ambiguous or not.
     */
    boolean isAmbiguousState(int state);

    /**
     * true if the class does not need any further initialisation.
     */
    boolean isStandard();

    /**
     * Gets data type description.
     */
    String getTypeDescription();

    /**
     * Gets character corresponding to a given state.
     */
    char getChar(int state);

    /**
     * Get a string code corresponding to a given state.
     */
    String getCode(int state);

    public abstract class Base implements DataType {

        protected int stateCount;
        protected String codeMap;
        protected int codeLength;
        protected int[][] mapCodeToStateSet;

        public String getCodeMap() {
            return codeMap;
        }

        @Override
        public int getStateCount() {
            return stateCount;
        }

        /**
         * implementation for single character per state encoding *
         */
        @Override
        public List<Integer> stringToState(String data) {
            List<Integer> sequence = new ArrayList<>();
            // remove spaces
            data = data.replaceAll("\\s", "");
            data = data.toUpperCase();
            if (codeMap == null) {
                if (data.contains(",")) {
                    // assume it is a comma separated string of integers
                    String[] strs = data.split(",");
                    for (String str : strs) {
                        try {
                            sequence.add(Integer.parseInt(str));
                        } catch (NumberFormatException e) {
                            sequence.add(-1);
                        }
                    }
                } else {
                    // assume it is a string where each character is a state
                    for (byte c : data.getBytes()) {
                        switch (c) {
                            case GAP_CHAR:
                            case MISSING_CHAR:
                                sequence.add(-1);
                                break;
                            default:
                                sequence.add(Integer.parseInt((char) c + ""));
                        }
                    }
                }
            } else {
                if (codeLength == 1) {
                    // single character codes
                    for (int i = 0; i < data.length(); i++) {
                        char cCode = data.charAt(i);
                        int stateCount = codeMap.indexOf(cCode);
                        if (stateCount < 0) {
                            throw new IllegalArgumentException("Unknown code found in sequence: " + cCode);
                        }
                        sequence.add(stateCount);
                    }
                } else if (codeLength > 1) {
                    // multi-character codes of fixed length
                    // use code map to resolve state codes
                    Map<String, Integer> map = new HashMap<>();
                    // fixed length code
                    for (int i = 0; i < codeMap.length(); i += codeLength) {
                        String code = codeMap.substring(i, i + codeLength);
                        map.put(code, i / codeLength);
                    }

                    for (int i = 0; i < data.length(); i += codeLength) {
                        String code = data.substring(i, i + codeLength).toUpperCase();
                        if (map.containsKey(code)) {
                            sequence.add(map.get(code));
                        } else {
                            throw new IllegalArgumentException("Unknown code found in sequence: " + code);
                        }
                    }
                } else {
                    // variable length code of strings
                    String[] codes = codeMap.toUpperCase().split(",");
                    for (String code : data.split(",")) {
                        boolean isFound = false;
                        for (int codeIndex = 0; codeIndex < codes.length; codeIndex++) {
                            if (code.equals(codes[codeIndex])) {
                                sequence.add(codeIndex);
                                isFound = true;
                                break;
                            }
                        }
                        if (!isFound) {
                            throw new RuntimeException("Could not find code " + code + " in codemap");
                        }
                    }
                }
            }
            return sequence;
        }

        @Override
        public String stateToString(List<Integer> nrOfStates) {
            int[] nrOfStates2 = new int[nrOfStates.size()];
            for (int i = 0; i < nrOfStates2.length; i++) {
                nrOfStates2[i] = nrOfStates.get(i);
            }
            return stateToString(nrOfStates2);
        }

        @Override
        public String stateToString(int[] nrOfStates) {
            StringBuffer buf = new StringBuffer();
            if (codeMap != null) {
                for (int state : nrOfStates) {
                    String code = codeMap.substring(state * codeLength, state * codeLength + codeLength);
                    buf.append(code);
                }
            } else {
                // produce a comma separated string of integers
                for (int i = 0; i < nrOfStates.length - 1; i++) {
                    buf.append(nrOfStates[i] + ",");
                }
                buf.append(nrOfStates[nrOfStates.length - 1] + "");
            }
            return buf.toString();
        }


        @Override
        public int[] getStatesForCode(int state) {
            return mapCodeToStateSet[state];
        }

        @Override
        public boolean[] getStateSet(int state) {
            boolean[] stateSet = new boolean[stateCount];
            int[] stateNumbers = getStatesForCode(state);
            for (int i : stateNumbers) {
                stateSet[i] = true;
            }
            return stateSet;
        }

        /**
         * Definition of ambiguousity
         * Non-ambiguous states: 0 <= X < stateCount
         * Ambiguous states: X >= stateCount
         * Missing data: X < 0
         */
        @Override
        public boolean isAmbiguousState(int state) {
            return (state < 0 || state >= stateCount);
        }

        @Override
        public boolean isStandard() {
            return true;
        }

        @Override
        public char getChar(int state) {
            return (char) (state + 'A');
        }

        @Override
        public String getCode(int state) {
            return String.valueOf(getChar(state));
        }

        @Override
        public String toString() {
            return getTypeDescription();
        }

        public Integer char2state(String character) {
            return stringToState(character).get(0);
        }
    }

}