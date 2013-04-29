package edu.rice.cs.bioinfo.programs.soranus.controllers;

import edu.rice.cs.bioinfo.library.epidemiology.transmissionMap.Snitkin2012.SnitkinEdge;
import edu.rice.cs.bioinfo.library.epidemiology.transmissionMap.TransMapResult;
import edu.rice.cs.bioinfo.library.graph.algorithms.spanningTree.GraphDisconnectedException;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.soranus.models.analysis.*;
import edu.rice.cs.bioinfo.programs.soranus.models.data.FirstPositiveDataProvider;
import edu.rice.cs.bioinfo.programs.soranus.models.data.Sequencing;
import edu.rice.cs.bioinfo.programs.soranus.models.data.SequencingsDataProvider;
import edu.rice.cs.bioinfo.programs.soranus.models.data.TraceDataProvider;
import edu.rice.cs.bioinfo.programs.soranus.viewModels.*;
import org.joda.time.LocalDate;

import java.util.*;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 4/1/13
 * Time: 3:15 PM
 * To change this template use File | Settings | File Templates.
 */
public class AnalysisController<SDR,TDR,FPDR, E, EX extends Throwable>
{


    public class AnalysisRecord extends edu.rice.cs.bioinfo.programs.soranus.viewModels.AnalysisRecord
    {
        public final DocumentVM DocumentVM;

        public AnalysisRecord(String title, DocumentVM vm) {
            super(title);
            DocumentVM = vm;
        }
    }

    private final WorkspaceVM<?,AnalysisController.AnalysisRecord> _workspaceVM;

    private final SequencingsDataProvider<SDR, E, EX> _sequencingsDataProvider;

    private final TraceDataProvider<TDR, E, EX> _traceDataProvider;

    private final FirstPositiveDataProvider<FPDR, E, EX> _firstPositiveDataProvider;

    public AnalysisController(WorkspaceVM<?,AnalysisController.AnalysisRecord> vm, SequencingsDataProvider<SDR, E, EX> sequencingsDataProvider, TraceDataProvider<TDR, E, EX> traceDataProvider,
                              FirstPositiveDataProvider<FPDR, E, EX> firstPositiveDataProvider) {
        _workspaceVM = vm;
        _sequencingsDataProvider = sequencingsDataProvider;
        _traceDataProvider = traceDataProvider;
        _firstPositiveDataProvider = firstPositiveDataProvider;



    }

    public void analysisRecordSelected(AnalysisRecord record)
    {
        _workspaceVM.setFocusDocument(record.DocumentVM);
    }

    public void performRecombDetectionUnderInfSite(SDR sequencingsDataRecord) throws EX
    {
        Set<Sequencing> sequencings = _sequencingsDataProvider.getSequencingToPatientMap(sequencingsDataRecord).keySet();

        LinkedList<Integer> sites = new LinkedList<Integer>();
        final String someSequence = sequencings.iterator().next().Sequence;

        for(int i = 0; i<someSequence.length(); i++)
        {
            sites.add(i);
        }


        RecombUnderInfiniteSites<Sequencing,Integer>.RecombProof proof = new RecombUnderInfiniteSites<Sequencing,Integer>()
        {
            @Override
            protected boolean isReferenceNucleotide(Sequencing sequence, Integer site)
            {
                return sequence.Sequence.charAt(site) == someSequence.charAt(0);
            }
        }.tryFindRecombinationProof(sequencings, sites);

        RecombResultVM vm = null;
        if(proof != null)
        {
            vm = new RecombResultVM(proof.Sequence1.toString(), proof.Sequence2.toString(), proof.Sequence3.toString(),
                                    proof.Sequence4.toString(),
                                    proof.Site1.toString(), proof.Site2.toString(),
                                    proof.Sequence1.Sequence.charAt(proof.Site1), proof.Sequence2.Sequence.charAt(proof.Site1),
                                    proof.Sequence3.Sequence.charAt(proof.Site1), proof.Sequence4.Sequence.charAt(proof.Site1),
                                    proof.Sequence1.Sequence.charAt(proof.Site2), proof.Sequence2.Sequence.charAt(proof.Site2),
                                    proof.Sequence3.Sequence.charAt(proof.Site2), proof.Sequence4.Sequence.charAt(proof.Site2));
        }

        AnalysisRecord ar = new AnalysisRecord("Recomb Result", vm);
        _workspaceVM.addAnalysis(ar);
        _workspaceVM.setFocusDocument(vm);



    }

    public void performSnitkinTransMapAnalysis(SDR sequencingsDataRecord, TDR traceDataRecord,
                                               FPDR firstPositiveDataRecord) throws EX {

        final Map<Sequencing, E> sequencings = _sequencingsDataProvider.getSequencingToPatientMap(sequencingsDataRecord);
        final Map<E, Map<LocalDate, Object>> patientTraces = _traceDataProvider.getTraces(traceDataRecord);
        final Map<E,LocalDate> firstPositiveDates = _firstPositiveDataProvider.getFirstPositiveDates(firstPositiveDataRecord);

        final TransMapResult<SnitkinEdge<E,Double>> transResult =
                new SnitkinTransMapInferrer<E>(patientTraces, firstPositiveDates, sequencings).inferMaps();
        final Set<SnitkinEdge<E,Double>> transMap = transResult.getSolutions().iterator().next();

        TransMapVM<E, SnitkinEdge<E,Double>> transMapViewModel =
                new TransMapVM<E, SnitkinEdge<E, Double>>(transMap) {
                    @Override
                    public E getSource(SnitkinEdge<E, Double> edge) {
                        return edge.getSource();
                    }

                    @Override
                    public E getDestination(SnitkinEdge<E, Double> edge) {
                        return edge.getDestination();
                    }

                    @Override
                    public String getNodeLabel(E source) {
                        return source.toString();
                    }

                    @Override
                    public String getEdgeLabel(SnitkinEdge<E,Double> edge)
                    {
                        return null;
                    }
                };

        AnalysisRecord ar = new AnalysisRecord("Snitkin Trans Map", transMapViewModel);
        _workspaceVM.addAnalysis(ar);
        _workspaceVM.setFocusDocument(transMapViewModel);


    }

    public void performMinSpanTreeAnalysisSnpMax2(SDR sequencingsDataRecord) throws EX, GraphDisconnectedException
    {
        final Map<Sequencing, E> sequeincingToSource = _sequencingsDataProvider.getSequencingToPatientMap(sequencingsDataRecord);

        final Map<E,String> sourceToDispalyString = new HashMap<E, String>();



        E prev = null;
        String greatestCommonPrefix = null;
        for(E entity : sequeincingToSource.values())
        {
            String entityString = entity.toString();
            sourceToDispalyString.put(entity, entityString);

            if(prev == null)
                greatestCommonPrefix  = entityString;
            else
                greatestCommonPrefix = greatestCommonPrefix(greatestCommonPrefix, entityString);

            prev = entity;
        }

        for(E entity : sourceToDispalyString.keySet())
        {
            sourceToDispalyString.put(entity,
                    sourceToDispalyString.get(entity).substring(greatestCommonPrefix.length()));
        }

        Iterable<Set<SequencingEdge<Sequencing>>> msts = new MinSpanTreesSnpInferrerMaxTwo<Sequencing>()
        {
            @Override
            protected Long getSnpDistance(Sequencing sequencing1, Sequencing sequencing2)
            {
                return new Long(sequencing1.getGeneticDistance(sequencing2));
            }
        }.inferMinTrees(sequeincingToSource.keySet());

                /*
        new MinSpanTreesSnpInferrerCombAsc<Sequencing>()
        {

            @Override
            protected Long getSnpDistance(Sequencing sequencing1, Sequencing sequencing2)
            {
                return  new Long(sequencing1.getGeneticDistance(sequencing2));
            }

        }.inferMinTrees(new HashSet(sequeincingToSource.keySet())); */

          /*      new MinSpanTreesSnpInferrerBacktrack<Sequencing>()
                {
                    @Override
                    protected Long getSnpDistance(Sequencing sequencing1, Sequencing sequencing2)
                    {
                        return new Long(sequencing1.getGeneticDistance(sequencing2));
                    }
                }.inferMinTrees(new HashSet(sequeincingToSource.keySet())); */

     /*   Set<SequencingEdge<Sequencing>> completeGraph = new CompleteGraphFactory<Sequencing,SequencingEdge<Sequencing>>()
        {

            @Override
            public SequencingEdge<Sequencing> makeEdge(Sequencing node1, Sequencing node2)
            {
                return new SequencingEdge<Sequencing>(node1, node2, (long)node1.getGeneticDistance(node2));
            }
        }.makeCompleteGraph(sequeincingToSource.keySet());

        completeGraph = new HashSet<SequencingEdge<Sequencing>>(IterableHelp.toList(IterableHelp.filter(completeGraph, new Predicate1<SequencingEdge<Sequencing>>()
        {
            public boolean execute(SequencingEdge<Sequencing> input)
            {
                return input.SnpDistance == 0;  //To change body of implemented methods use File | Settings | File Templates.
            }
        })));

        TreeVM<Sequencing,SequencingEdge<Sequencing>> vmm = new
                TreeVM<Sequencing,SequencingEdge<Sequencing>>(completeGraph)
                {

                    @Override
                    public Tuple<Sequencing, Sequencing> getNodesOfEdge(SequencingEdge<Sequencing> edge)
                    {
                        return new Tuple<Sequencing, Sequencing>(edge.Sequencing1, edge.Sequencing2);
                    }

                    @Override
                    public String getNodeLabel(Sequencing sequencing)
                    {
                        return sourceToDispalyString.get(sequeincingToSource.get(sequencing));
                    }

                    @Override
                    public String getEdgeLabel(SequencingEdge<Sequencing> edge)
                    {
                        return edge.SnpDistance + "";
                    }
                };

        AnalysisRecord rr = new AnalysisRecord("SNP Graph", vmm);
        _workspaceVM.addAnalysis(rr);
        _workspaceVM.setFocusDocument(vmm);      */

        final Set<SequencingEdge<Sequencing>> commonEdges = new HashSet<SequencingEdge<Sequencing>>();
        final Set<SequencingEdge<Sequencing>> firstEdgesOnly = new HashSet<SequencingEdge<Sequencing>>();
        Iterator<Set<SequencingEdge<Sequencing>>> mstElements = msts.iterator();
        Set<SequencingEdge<Sequencing>> mst1 = mstElements.next();
        Set<SequencingEdge<Sequencing>> mst2 = mstElements.next();

        for(SequencingEdge<Sequencing> edge1 : mst1)
        {
            for(SequencingEdge<Sequencing> edge2 : mst2)
            {
                if(sameEdge(edge1, edge2))
                {
                    commonEdges.add(edge1);
                    commonEdges.add(edge2);
                    break;
                }
            }

            if(!commonEdges.contains(edge1))
                firstEdgesOnly.add(edge1);
        }

        Set<SequencingEdge<Sequencing>> secondEdgesOnly = new HashSet<SequencingEdge<Sequencing>>(mst2);
        secondEdgesOnly.removeAll(commonEdges);


        Set<SequencingEdge<Sequencing>> allEdges = new HashSet<SequencingEdge<Sequencing>>(commonEdges);
        allEdges.addAll(firstEdgesOnly);
        allEdges.addAll(secondEdgesOnly);
        TreeVM<Sequencing,SequencingEdge<Sequencing>> vm = new
                TreeVM<Sequencing,SequencingEdge<Sequencing>>(allEdges)
                {
                    @Override
                    public int getEdgeColor(SequencingEdge<Sequencing> edge)
                    {
                        if(commonEdges.contains(edge))
                            return 0;

                        if(firstEdgesOnly.contains(edge))
                            return 1;

                        return 2;
                    }

                    @Override
                    public Tuple<Sequencing, Sequencing> getNodesOfEdge(SequencingEdge<Sequencing> edge)
                    {
                        return new Tuple<Sequencing, Sequencing>(edge.Sequencing1, edge.Sequencing2);
                    }

                    @Override
                    public String getNodeLabel(Sequencing sequencing)
                    {
                        return sourceToDispalyString.get(sequeincingToSource.get(sequencing));
                    }

                    @Override
                    public String getEdgeLabel(SequencingEdge<Sequencing> edge)
                    {
                        return edge.SnpDistance + "";
                    }
                };


        AnalysisRecord ar = new AnalysisRecord("Min Span Tree (SNP)", vm);
        _workspaceVM.addAnalysis(ar);
        _workspaceVM.setFocusDocument(vm);
    }

    private boolean sameEdge(SequencingEdge<Sequencing> edge1, SequencingEdge<Sequencing> edge2)
    {
        return (edge1.Sequencing1 == edge2.Sequencing1 && edge1.Sequencing2 == edge2.Sequencing2) ||
                (edge1.Sequencing1 == edge2.Sequencing2 && edge1.Sequencing2 == edge2.Sequencing1);
    }

    public void performNeighborJoiningAnalysis(SDR sequencingsDataRecord) throws EX {


        final Map<Sequencing, E> sequeincingToSource = _sequencingsDataProvider.getSequencingToPatientMap(sequencingsDataRecord);


        NeighborJoining nj = new NeighborJoining();
        NeighborJoining.Graph joinTree = nj.performJoin(sequeincingToSource.keySet());




        NeighborJoiningVM<Sequencing,NeighborJoining.Edge> joinVM = new NeighborJoiningVM<Sequencing,NeighborJoining.Edge>(joinTree.Nodes, joinTree.Edges)
        {
            @Override
            public Tuple<Sequencing, Sequencing> getNodesOfEdge(NeighborJoining.Edge edge) {
                return new Tuple<Sequencing, Sequencing>(edge.Node1, edge.Node2);
            }

            @Override
            public String getNodeLabel(Sequencing node) {
                if(sequeincingToSource.containsKey(node))
                {
                    return sequeincingToSource.get(node).toString();
                }
                else
                {
                    return "";
                }
            }

            @Override
            public String getEdgeLabel(NeighborJoining.Edge edge) {
                return edge.Distance.toString();
            }
        };

        AnalysisRecord ar = new AnalysisRecord("Neighbor Joining", joinVM);
        _workspaceVM.addAnalysis(ar);
        _workspaceVM.setFocusDocument(joinVM);


    }

    private String greatestCommonPrefix(String a, String b) {
        int minLength = Math.min(a.length(), b.length());
        for (int i = 0; i < minLength; i++) {
            if (a.charAt(i) != b.charAt(i)) {
                return a.substring(0, i);
            }
        }
        return a.substring(0, minLength);
    }
}
