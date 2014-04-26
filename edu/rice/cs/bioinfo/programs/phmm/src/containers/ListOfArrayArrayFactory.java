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

//package containers;
//
//import java.util.ArrayList;
//
//public class ListOfArrayArrayFactory<Type> implements IArrayFactory<IListOfArrays<Type>> {
//	@Override
//	public IArray<IListOfArrays<Type>> make(final int length) {
//		return new IArray<IListOfArrays<Type>>() {
//			ArrayList<IListOfArrays<Type>> list = new ArrayList<IListOfArrays<Type>>(length);
//
//			@Override
//			public IListOfArrays<Type> get(int index) {
//				return list.get(index);
//			}
//
//			@Override
//			public void set(int index, IListOfArrays<Type> type) {
//				list.set(index, type);
//			}
//
//			@Override
//			public int getLength() {
//				return length;
//			}
//		};
//	}
//}
