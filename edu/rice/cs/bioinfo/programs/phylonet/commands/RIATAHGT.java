package edu.rice.cs.bioinfo.programs.phylonet.commands;

import com.sun.org.apache.xalan.internal.xsltc.compiler.Template;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.*;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.ast.*;
import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.Proc3;
import edu.rice.cs.bioinfo.programs.phylonet.algos.riatahgt.EventBootstrap;
import edu.rice.cs.bioinfo.programs.phylonet.algos.riatahgt.HgtEvent;
import edu.rice.cs.bioinfo.programs.phylonet.algos.riatahgt.HgtScenario;
import edu.rice.cs.bioinfo.programs.phylonet.algos.riatahgt.RiataHgt;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.io.HgtReader;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io.NewickReader;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import javax.xml.transform.Result;
import java.io.*;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 10/6/11
 * Time: 1:28 PM
 * To change this template use File | Settings | File Templates.
 */
public class RIATAHGT extends CommandBaseFileOut
{
    NetworkNonEmpty _speciesTree;

    LinkedList<NetworkNonEmpty> _geneTrees = new LinkedList<NetworkNonEmpty>();

    private Parameter _expandedOutputParam = null;

    private String _nodePrefixValue = null;

    private boolean _refined = true;

    private boolean _collapsed = true;

    RIATAHGT(SyntaxCommand motivatingCommand, ArrayList<Parameter> params, Map<String, NetworkNonEmpty> sourceIdentToNetwork, Proc3<String, Integer, Integer> errorDetected) {
        super(motivatingCommand, params, sourceIdentToNetwork, errorDetected);
    }

    protected int getMinNumParams()
    {
        return 2;
    }

    protected int getMaxNumParams()
    {
        return 6;
    }

    public boolean checkParamsForCommand() {

         boolean noError = true;

         int expectedSpeciesTreeParamIndex = 0;

         for(int i = 0; i<=3 && i<this.params.size(); i++)
         {
             Parameter ithParam = params.get(i);

             String ithParamValue = ithParam.execute(GetSimpleParamValue.Singleton, null);

             if(ithParamValue != null)
             {
                 if(ithParamValue.toLowerCase().equals("-p"))
                 {
                     if(_nodePrefixValue == null)
                     {
                         if(i + 1 < this.params.size())
                         {
                             Parameter ithPlusOneParam = this.params.get(i + 1);
                             _nodePrefixValue = ithPlusOneParam.execute(GetSimpleParamValue.Singleton, null);
                             expectedSpeciesTreeParamIndex = i + 2;

                             if(_nodePrefixValue == null)
                             {
                                 noError = false;
                                 errorDetected.execute("Invalid node prefix.", ithPlusOneParam.getLine(), ithPlusOneParam.getColumn());
                             }
                         }
                         else
                         {
                             noError = false;
                             errorDetected.execute("Expected value after switch -p.", ithParam.getLine(), ithParam.getColumn());
                         }
                     }
                     else
                     {
                         noError = false;
                         errorDetected.execute("Duplicate switch -p", ithParam.getLine(), ithParam.getColumn());
                     }
                 }
                 else if(ithParamValue.toLowerCase().equals("-e"))
                 {
                     _expandedOutputParam = this.params.get(i);
                     expectedSpeciesTreeParamIndex = i + 1;
                 }
                 else if(ithParamValue.toLowerCase().equals("-u"))
                 {
                     _refined = false;
                     _collapsed = false;
                     expectedSpeciesTreeParamIndex = i + 1;
                 }
             }
         }


         ParameterIdentSet geneTreesParam = null;
         ParameterIdent speciesTreeParam = null;



         noError = this.params.get(expectedSpeciesTreeParamIndex).execute(new ParameterAlgo<Boolean, Boolean, RuntimeException>() {
             public Boolean forIdentifier(ParameterIdent parameter, Boolean o) throws RuntimeException {
                 return o;  //To change body of implemented methods use File | Settings | File Templates.
             }

             public Boolean forIdentList(ParameterIdentList parameterIdentList, Boolean aBoolean) throws RuntimeException {
                 return null;  //To change body of implemented methods use File | Settings | File Templates.
             }

             public Boolean forQuote(ParameterQuote parameter, Boolean o) throws RuntimeException {
                 return o;  //To change body of implemented methods use File | Settings | File Templates.
             }

             public Boolean forTaxonSetList(ParameterTaxonSetList parameterTaxonSetList, Boolean o) throws RuntimeException {
                 return error(parameterTaxonSetList);
             }

             public Boolean forIdentSet(ParameterIdentSet parameterIdentSet, Boolean o) throws RuntimeException {
                return error(parameterIdentSet);
             }

             public Boolean forTaxaMap(ParameterTaxaMap parameterTaxaMap, Boolean aBoolean) throws RuntimeException {
                 return error(parameterTaxaMap);
             }

             private Boolean error(Parameter p)
             {
                errorDetected.execute("Expected species tree identifier.", p.getLine(), p.getColumn());
                return false;  //To change body of implemented methods use File | Settings | File Templates.
             }

         }, noError);

         if(noError)
         {
             speciesTreeParam = this.assertParameterIdent(params.get(expectedSpeciesTreeParamIndex));

             if(speciesTreeParam == null)
             {
                 noError = false;
             }

         }

         noError = this.params.get(expectedSpeciesTreeParamIndex + 1).execute(new ParameterAlgo<Boolean, Boolean, RuntimeException>()
         {
             public Boolean forIdentifier(ParameterIdent parameter, Boolean aBoolean) throws RuntimeException {
                 error(parameter);
                 return false;  //To change body of implemented methods use File | Settings | File Templates.
             }

             public Boolean forIdentList(ParameterIdentList parameterIdentList, Boolean aBoolean) throws RuntimeException {
                 error(parameterIdentList);
                 return false;  //To change body of implemented methods use File | Settings | File Templates.
             }

             public Boolean forQuote(ParameterQuote parameter, Boolean aBoolean) throws RuntimeException {
                 error(parameter);
                 return false;  //To change body of implemented methods use File | Settings | File Templates.
             }

             public Boolean forTaxonSetList(ParameterTaxonSetList parameterTaxonSetList, Boolean aBoolean) throws RuntimeException {
                 error(parameterTaxonSetList);
                 return false;  //To change body of implemented methods use File | Settings | File Templates.
             }

             public Boolean forIdentSet(ParameterIdentSet parameterIdentSet, Boolean aBoolean) throws RuntimeException {
                 return aBoolean;  //To change body of implemented methods use File | Settings | File Templates.
             }

             public Boolean forTaxaMap(ParameterTaxaMap parameterTaxaMap, Boolean aBoolean) throws RuntimeException {
                 error(parameterTaxaMap);
                 return false;  //To change body of implemented methods use File | Settings | File Templates.
             }

             private void error(Parameter parameter) {
                 errorDetected.execute("Expected gene trees identifier set. (e.g. {tree1, tree2}).", parameter.getLine(), parameter.getColumn());
             }

         }, noError);

         if(noError)
         {
             geneTreesParam = (ParameterIdentSet) this.params.get(expectedSpeciesTreeParamIndex + 1);
         }


         if(noError)
         {
             return checkContext(sourceIdentToNetwork, errorDetected, geneTreesParam, speciesTreeParam);
         }
         else
         {
             return false;
         }
    }

    private boolean checkContext(Map<String, NetworkNonEmpty> sourceIdentToNetwork, final Proc3<String, Integer, Integer> errorDetected,
                                 ParameterIdentSet geneTrees, ParameterIdent speciesTree)
    {
        boolean noError = true;

        noError = this.assertNetworkExists(speciesTree);

        if(noError)
        {
            _speciesTree = sourceIdentToNetwork.get(speciesTree.execute(GetSimpleParamValue.Singleton, null));
        }

        for(String geneTreeIdent : geneTrees.Elements)
        {
            noError = this.assertNetworkExists(geneTreeIdent, geneTrees.getLine(), geneTrees.getColumn());
            if(noError)
            {
                _geneTrees.add(sourceIdentToNetwork.get(geneTreeIdent));
            }
        }


        return noError;
    }

    @Override
    protected String produceResult() {

        final StringBuffer result = new StringBuffer();
        STITree<Double> speciesTree = new STITree<Double>(true);
		List<STITree<Double>> geneTrees = new LinkedList<STITree<Double>>();

		String prefix = "";
		boolean show_compact = true;

		       /*
		BufferedReader br;
		if (infile == null) {
			System.out.println("Enter a species tree, and gene trees on separate lines:");
			br = new BufferedReader(new InputStreamReader(System.in));
		}
		else {
			try {
				br = new BufferedReader(new FileReader(infile));
			}
			catch (FileNotFoundException e) {
				System.err.println(e.getMessage());
				return;
			}
		} */

		// Read the input for species and gene trees.
        SingleLinePrinter printer = new SingleLinePrinter();
        printer.setSupportTransformer(new TransformSupportToBase100());
		String firstLine = printer.toString(_speciesTree);
		List<String> nextLines = new LinkedList<String>();


			try {
				NewickReader nr = new NewickReader(new StringReader(firstLine));


				// speciesTree = (STITree<Object>) nr.readTree();
				nr.readTree(speciesTree);
				Trees.removeBinaryNodes(speciesTree);

				for (NetworkNonEmpty g : _geneTrees) {
                    String s = printer.toString(g);
					nr = new NewickReader(new StringReader(s));

					// STITree<Object> gt = (STITree<Object>) nr.readTree();
					STITree<Double> gt = new STITree<Double>(true);
					nr.readTree(gt);

					Trees.removeBinaryNodes(gt);
					geneTrees.add(gt);
				}
			}
			catch(Exception e)
            {
                throw new RuntimeException(e);
            }


		// Compute HGT events and print to the output.
		int gt_idx = 0;
		List<List<HgtEvent>> allGenesEvents	= new LinkedList<List<HgtEvent>>();	// Store all events detected by RIATA to build the consensus network.
		STITree<Double> copy = null;

		for (STITree<Double> gt : geneTrees) {
			RiataHgt hgt = new RiataHgt();

            if(_nodePrefixValue != null)
            {
                hgt.setPrefix(_nodePrefixValue);
            }

			if (!_collapsed) {
				hgt.disableCollapse();
			}
			if (!_refined) {
				hgt.disableRefine();
			}

			copy = new STITree<Double>(speciesTree.getRoot());

			// Keep only leaves that appear in both trees.
			Trees.pruneLeaves(copy, gt);
			hgt.computeHgt(copy, gt);

			// Store all events found for this gene tree for computing the consensus network.
			List<HgtEvent> events  = new LinkedList<HgtEvent>();

			for (STINode<List<HgtScenario>> node : hgt.getSolutionTree().getNodes()) {
				for (HgtScenario hs : node.getData()) {
					for (HgtEvent he : hs.getEvents()) {
						if (!events.contains(he)) {
							events.add(he);
						}
					}
				}
			}

			allGenesEvents.add(events);

			// Compute event bootstrap.
			EventBootstrap eb = new EventBootstrap(copy, gt, hgt.getSolutionTree());
			eb.computeBootstrap();

			if (_expandedOutputParam == null) {
				result.append("species tree: " + copy.toString()+ "\n");
				result.append("gene tree: " + gt.toString()+ "\n");

				// Add bootstrap values to events.

				StringBuffer summary = new StringBuffer(hgt.toString());
				for (int i = 0; i < eb.getBootstraps().size(); i++) {
					if (eb.getBootstraps().get(i) == TNode.NO_DISTANCE) {
						continue;	// Skip events with no bootstraps.
					}

					String delim = eb.getEvents().get(i).toString();
					String val = "(" + eb.getBootstraps().get(i) + ")";

					int from = 0;
					while ((from = summary.indexOf(delim, from)) != -1) {
						from += delim.length();
						summary.insert(from, val);
					}
				}

				result.append(summary.toString()+ "\n");	// Print summary of events and solutions.
				result.append("*****************************************************************************************************"+ "\n");



				gt_idx++;
			}
			else {
				// Display in the non-compact format.
				result.append("species tree: " + copy.toString()+ "\n");
				result.append("gene tree: " + gt.toString()+ "\n");

				List<HgtScenario> scenarios = hgt.enumerateSolutions();
				int count = 1;
				for (HgtScenario hs : scenarios) {
					result.append("Solution " + count+ "\n");

					// Display in the eNewick format.
					// ps.print(hs.toString());
					String temp = new String();
					temp = copy.toString() + "\n";
					for (HgtEvent event : hs.getEvents()) {
						temp += event.getSourceEdge().getName() + " -> " + event.getDestEdge().getName() + "\n";
					}

					try {
						HgtReader<Object> hr = new HgtReader<Object>(new StringReader(temp));
						edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network <Object> net = hr.readHgt();
                        Network newNet = NetworkTransformer.fromClassicNetwork(net);

                        newNet.execute(new NetworkAlgo<Object, Object, RuntimeException>() {
                            public Object forNetworkEmpty(NetworkEmpty networkEmpty, Object o) throws RuntimeException {

                                return null;
                            }

                            public Object forNetworkNonEmpty(NetworkNonEmpty networkNonEmpty, Object o) throws RuntimeException {
                                result.append(new SingleLinePrinter().toString(networkNonEmpty));
                                return null;
                            }
                        }, null);
                        result.append("\n");


				//		ExNewickWriter<Object> netWriter = new ExNewickWriter<Object>(new PrintWriter(ps));
			//			netWriter.writeNetwork(net);
			//			ps.println();
					}
					catch (Exception e) {
						throw new RuntimeException(e);
					}

					count++;
					if (count <= scenarios.size()) {
						result.append("-----------------------------------------------------------------------------------------------------\n");
					}
				}

				result.append("*****************************************************************************************************\n");
			}
		}

		// Find the consensus network.
		List<HgtEvent> consensusEvents = RiataHgt.getConsensusNetwork(allGenesEvents);

		result.append("Consensus network for this set of gene trees\n");
		result.append(copy.toString() + "\n");
		for (HgtEvent he : consensusEvents) {
			result.append(he.getSourceEdge().getName() + " -> " + he.getDestEdge().getName() + "\n");
		}




        return result.toString();
    }
}
