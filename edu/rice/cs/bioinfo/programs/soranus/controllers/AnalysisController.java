package edu.rice.cs.bioinfo.programs.soranus.controllers;

import edu.rice.cs.bioinfo.library.epidemiology.transmissionMap.Snitkin2012.SnitkinEdge;
import edu.rice.cs.bioinfo.library.epidemiology.transmissionMap.TransMapResult;
import edu.rice.cs.bioinfo.library.graph.algorithms.spanningTree.GraphDisconnectedException;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.soranus.models.analysis.MinSpanTreesSnpInferrerCombAsc;
import edu.rice.cs.bioinfo.programs.soranus.models.analysis.NeighborJoining;
import edu.rice.cs.bioinfo.programs.soranus.models.analysis.SequencingEdge;
import edu.rice.cs.bioinfo.programs.soranus.models.analysis.SnitkinTransMapInferrer;
import edu.rice.cs.bioinfo.programs.soranus.models.data.FirstPositiveDataProvider;
import edu.rice.cs.bioinfo.programs.soranus.models.data.Sequencing;
import edu.rice.cs.bioinfo.programs.soranus.models.data.SequencingsDataProvider;
import edu.rice.cs.bioinfo.programs.soranus.models.data.TraceDataProvider;
import edu.rice.cs.bioinfo.programs.soranus.viewModels.*;
import org.joda.time.LocalDate;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

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

    public void performMinSpanTreeAnalysisSnp(SDR sequencingsDataRecord) throws EX, GraphDisconnectedException
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

        Set<Set<SequencingEdge<Sequencing>>> msts =

        new MinSpanTreesSnpInferrerCombAsc<Sequencing>()
        {

            @Override
            protected Long getSnpDistance(Sequencing sequencing1, Sequencing sequencing2)
            {
                return  new Long(sequencing1.getGeneticDistance(sequencing2));
            }

        }.inferMinTrees(new HashSet(sequeincingToSource.keySet()));

          /*      new MinSpanTreesSnpInferrerBacktrack<Sequencing>()
                {
                    @Override
                    protected Long getSnpDistance(Sequencing sequencing1, Sequencing sequencing2)
                    {
                        return new Long(sequencing1.getGeneticDistance(sequencing2));
                    }
                }.inferMinTrees(new HashSet(sequeincingToSource.keySet())); */

        for(Set<SequencingEdge<Sequencing>> mst : msts)
        {

            TreeVM<Sequencing,SequencingEdge<Sequencing>> vm = new
                    TreeVM<Sequencing,SequencingEdge<Sequencing>>(mst)
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
            AnalysisRecord ar = new AnalysisRecord("Min Span Tree (SNP)", vm);
            _workspaceVM.addAnalysis(ar);
            _workspaceVM.setFocusDocument(vm);
        }
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
