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

package edu.rice.cs.bioinfo.library.language.pyson._1_0.ast;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 9/6/11
 * Time: 2:55 PM
 * To change this template use File | Settings | File Templates.
 */
public interface BlockAlgo<R, T, E extends Exception> {

    public R forTreesBlock(TreesBlockBody treesBlock, T input) throws E;

    public R forNetworksBlock(NetworksBlockBody treesBlock, T input) throws E;

    public R forPhylonetBlockBody(PhyloNetBlockBody phyloBlock, T input) throws E;

    public R forDataBlock(DataBlockBody dataBlock, T input) throws E;

    public R forMorphDataBlock(MorphDataBlockBody morphDataBlock, T input) throws E;
}
