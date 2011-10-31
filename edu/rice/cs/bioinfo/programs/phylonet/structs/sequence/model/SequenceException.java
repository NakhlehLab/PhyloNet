package edu.rice.cs.bioinfo.programs.phylonet.structs.sequence.model;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 10/31/11
 * Time: 3:33 PM
 * To change this template use File | Settings | File Templates.
 */


    /**
     * An <code>Exception</code> type used by SequenceAlignment whenever an
     * invalidly-formatted sequence is parsed.
     */
    public class SequenceException extends Exception {

        public SequenceException(String message) {
            super(message);
        }

        public SequenceException(String message, Throwable cause) {
            super(message, cause);
        }

        public SequenceException(Throwable cause) {
            super(cause);
        }

    }


