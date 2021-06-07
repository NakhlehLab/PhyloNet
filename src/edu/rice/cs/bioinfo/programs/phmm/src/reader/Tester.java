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

//package reader;
//
//import java.util.ArrayList;
//
//import be.ac.ulg.montefiore.run.jahmm.Hmm;
//import be.ac.ulg.montefiore.run.jahmm.ObservationInteger;
//import be.ac.ulg.montefiore.run.jahmm.OpdfIntegerFactory;
//import be.ac.ulg.montefiore.run.jahmm.learn.BaumWelchLearner;
//
//public class Tester {
//
//	/**
//	 * @param args
//	 */
//	public static void main(String[] args) {
//
//		ArrayList<String> files = new ArrayList<String>();
//		files.add("myseq.txt");
//		files.add("myseq2.txt");
//		files.add("myseq3.txt");
//		files.add("myseq4.txt");
//		files.add("myseq5.txt");
//
//		try {
//			Parser a = new Parser("basicInfo.txt");
//			a.multiFileParser(files);
//
//
//
//		// Build HMM
//		Hmm<ObservationInteger > hmm = new Hmm<ObservationInteger>(5 , new OpdfIntegerFactory(256));
//
//		BaumWelchLearner bwl = new BaumWelchLearner();
//		bwl.learn(hmm, a.getObsList());
//
//		System.out.println("\n\n\n" + hmm);
//		} catch (Exception e) {
//			e.printStackTrace();
//		}
//
//
//	}
//
//}
