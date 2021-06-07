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


