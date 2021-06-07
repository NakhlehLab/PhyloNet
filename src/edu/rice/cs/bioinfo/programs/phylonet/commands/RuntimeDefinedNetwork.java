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

import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.NetworkNonEmpty;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 11/15/12
 * Time: 11:33 AM
 * To change this template use File | Settings | File Templates.
 */
public final class RuntimeDefinedNetwork extends NetworkNonEmpty
{
    public static final RuntimeDefinedNetwork Singleton = new RuntimeDefinedNetwork();

    private RuntimeDefinedNetwork() {
        super(null, null, null, null);

    }
}
