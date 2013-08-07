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
