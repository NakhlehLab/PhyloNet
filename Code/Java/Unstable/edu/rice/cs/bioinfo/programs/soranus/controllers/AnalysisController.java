package edu.rice.cs.bioinfo.programs.soranus.controllers;

import edu.rice.cs.bioinfo.library.epidemiology.transmissionMap.Snitkin2012.SnitkinEdge;
import edu.rice.cs.bioinfo.library.epidemiology.transmissionMap.TransMapResult;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.soranus.models.analysis.NeighborJoining;
import edu.rice.cs.bioinfo.programs.soranus.models.analysis.SnitkinTransMapInferrer;
import edu.rice.cs.bioinfo.programs.soranus.models.data.FirstPositiveDataProvider;
import edu.rice.cs.bioinfo.programs.soranus.models.data.Sequencing;
import edu.rice.cs.bioinfo.programs.soranus.models.data.SequencingsDataProvider;
import edu.rice.cs.bioinfo.programs.soranus.models.data.TraceDataProvider;
import edu.rice.cs.bioinfo.programs.soranus.viewModels.DocumentVM;
import edu.rice.cs.bioinfo.programs.soranus.viewModels.NeighborJoiningVM;
import edu.rice.cs.bioinfo.programs.soranus.viewModels.TransMapVM;
import edu.rice.cs.bioinfo.programs.soranus.viewModels.WorkspaceVM;
import org.joda.time.LocalDate;

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
                };

        AnalysisRecord ar = new AnalysisRecord("Snitkin Trans Map", transMapViewModel);
        _workspaceVM.addAnalysis(ar);
        _workspaceVM.setFocusDocument(transMapViewModel);


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
}
