package edu.rice.cs.bioinfo.library.epidemiology.transmissionMap.Snitkin2012;

import edu.rice.cs.bioinfo.library.epidemiology.transmissionMap.TransMapInferrer;
import edu.rice.cs.bioinfo.library.programming.*;
import edu.uci.ics.jung.graph.DirectedGraph;
import edu.uci.ics.jung.graph.DirectedSparseGraph;
import net.algowiki.edmonds._20121219.AdjacencyList;
import net.algowiki.edmonds._20121219.Edge;
import net.algowiki.edmonds._20121219.Edmonds;
import net.algowiki.edmonds._20121219.Node;
import org.joda.time.Days;
import org.joda.time.LocalDate;

import java.util.*;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 12/19/12
 * Time: 12:56 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class SnitkinTransMapInferrer<P, S,D extends Comparable<D>> implements TransMapInferrer<P,SnitkinEdge<P,D>>
{
    private final Set<P> _patients;

    private final List<S> _genomes;

    private final Func2<S,S,D> _getGeneticDistance;

    private final Map<P, Map<LocalDate,Object>> _patientTraces;

    private final Map<P, LocalDate> _patientToFirstPositiveDate;

    private final Map<S,P> _genomeToPatient;

    private P _rootPatient;

    private final Func2<D,D,D> _subtract = new Func2<D, D, D>() {
        public D execute(D a, D b) {
            return subtract(a, b);
        }
    };

    private DirectedGraph<P,SnitkinEdge<P,D>> makeEmptyGraph()
    {
        return new DirectedSparseGraph<P, SnitkinEdge<P, D>>();
    }

    public SnitkinTransMapInferrer(P rootPatient, List<S> genomes, Func2<S, S, D> getGeneticDistance,
                                   Map<P, Map<LocalDate,Object>> patientTraces,
                                   Map<P, LocalDate> patientToFirstPositiveDate,
                                   Map<S,P> genomeToPatient) {

        _patients = patientToFirstPositiveDate.keySet();

        if(!_patients.contains(rootPatient))
        {
            throw new IllegalArgumentException("root patient is not in the set of patients");
        }


        _genomes = genomes;
        _rootPatient = rootPatient;
        _getGeneticDistance = getGeneticDistance;
        _patientTraces = patientTraces;
        _patientToFirstPositiveDate = patientToFirstPositiveDate;
        _genomeToPatient = genomeToPatient;
    }


    protected abstract D subtract(D a, D b);

    protected abstract D makeZero();

    protected abstract D makeMax();

    protected abstract D makeDistance(int intDistance);

    protected abstract SnitkinEdge<P,D> makeEdge(P source, P destination, D geneticDistance, D epidemiologicalDistance);

    public Set<DirectedGraph<P, SnitkinEdge<P,D>>> inferMaps()
    {
        Map<P, Map<P, D>> geneticDistanceBetweenPatients = computeGeneticDistanceBetweenPatients();

        AdjacencyList digraph = new AdjacencyList();
        Map<Node,P> nodeToPatient = new HashMap<Node, P>();
        Map<P,Node> patientToNode = new HashMap<P, Node>();
        for(P patient : _patients)
        {
            Node pNode = new Node(nodeToPatient.size());
            nodeToPatient.put(pNode, patient);
            patientToNode.put(patient, pNode);
        }

        Set<SnitkinEdge<P,D>> completePatientGraphEdges = new HashSet<SnitkinEdge<P, D>>();
        for(P patient1 : _patients)
        {
            Node p1Node = patientToNode.get(patient1);
            for(P patient2 : _patients)
            {
                if(patient1 != patient2)
                {
                    Node p2Node = patientToNode.get(patient2);
                    D geneticDistance = geneticDistanceBetweenPatients.get(patient1).get(patient2);
                    D epidemiologicalDistance = getEpidemiologicalDistance(patient1,  patient2);
                    SnitkinEdge<P,D> edge = makeEdge(patient1, patient2, geneticDistance, epidemiologicalDistance);
                    completePatientGraphEdges.add(edge);
                    D distance = edge.getDistance();
                    digraph.addEdge(p1Node, p2Node, distance);
                    System.out.println(patient1 + " " + patient2 + " " + distance);

                }
            }
        }

        Node root = patientToNode.get(_rootPatient);

        Edmonds<D> edmonds = new Edmonds<D>();
        AdjacencyList<D> minimumSpan = edmonds.getMinBranching(root, digraph, _subtract, makeZero());

        DirectedGraph<P,SnitkinEdge<P,D>> resultGraph = makeEmptyGraph();
        for(P patient : _patients)
        {
            resultGraph.addVertex(patient);
        }

        for(Edge<D> minEdge : minimumSpan.getAllEdges())
        {
            P from = nodeToPatient.get(minEdge.from);
            P to = nodeToPatient.get(minEdge.to);
            System.out.println(from + " " + to);
            SnitkinEdge<P,D> resultEdge = findEdge(from, to, completePatientGraphEdges);
            resultGraph.addEdge(resultEdge, resultEdge.getSource(), resultEdge.getDestination());
        }

        return new HashSet<DirectedGraph<P, SnitkinEdge<P, D>>>(Arrays.asList(resultGraph));


    }

    private Map<P, Map<P, D>> computeGeneticDistanceBetweenPatients()
    {
        Map<P,Map<P,D>> geneticDistanceBetweenPatients = new HashMap<P, Map<P, D>>();
        for(P patient : _patients)
        {
            geneticDistanceBetweenPatients.put(patient, new HashMap<P, D>());
        }

        for(S genome1 : _genomes)
        {
            P patientOfGenome1 = _genomeToPatient.get(genome1);
            Map<P,D> pateint1Distances = geneticDistanceBetweenPatients.get(patientOfGenome1);

            for(S genome2 : _genomes)
            {
                P patientOfGenome2 = _genomeToPatient.get(genome2);

                if(patientOfGenome1.equals(patientOfGenome2))
                    continue;

                if(!_patients.contains(patientOfGenome1) || !_patients.contains(patientOfGenome2))
                    continue;

                D geneticDistance = _getGeneticDistance.execute(genome1, genome2);
                if(geneticDistance == null)
                    throw new RuntimeException("Genetic distance may not be null.");

                if(pateint1Distances.containsKey(patientOfGenome2))
                {
                    D currentMin = pateint1Distances.get(patientOfGenome2);
                    if(currentMin.compareTo(geneticDistance) > 0)
                    {
                        pateint1Distances.put(patientOfGenome2, geneticDistance);
                    }
                }
                else
                {
                    pateint1Distances.put(patientOfGenome2, geneticDistance);
                }

            }
        }
        return geneticDistanceBetweenPatients;
    }

    private D getEpidemiologicalDistance(P donor, P recipient)
    {
        LocalDate donorFirstPositive = _patientToFirstPositiveDate.get(donor);
        LocalDate recipientFirstPositive = _patientToFirstPositiveDate.get(recipient);

        LocalDate firstOverlap;
        try
        {
            firstOverlap = findFirstOverlap(donor, recipient);
        }
        catch (NoOverlapException e)
        {
            return makeMax();
        }

        if(recipientFirstPositive.isBefore(firstOverlap))
        {
            return makeMax();
        }
        else if(donorFirstPositive.isBefore(recipientFirstPositive))
        {
            LocalDate previousOverlap = findPreviousOverlap(donor, recipient, recipientFirstPositive);
            int numDaysBetween = Days.daysBetween(previousOverlap, recipientFirstPositive).getDays();
            return makeDistance(numDaysBetween);
        }
        else if(recipientFirstPositive.isBefore(donorFirstPositive))
        {
            int numDaysBetween = Days.daysBetween(recipientFirstPositive, donorFirstPositive).getDays();
            return makeDistance(numDaysBetween);
        }
        else
        {
            throw new IllegalArgumentException("Patient distance undefined");
        }
    }

    private LocalDate findPreviousOverlap(P patient1, P patient2, LocalDate onOrBefore) {

        Map<LocalDate,Object> patient1Locations = _patientTraces.get(patient1);
        Map<LocalDate,Object> patient2Locations = _patientTraces.get(patient2);

        List<LocalDate> patient1TraceDates = new ArrayList(patient1Locations.keySet());
        Collections.sort(patient1TraceDates);
        Collections.reverse(patient1TraceDates);
        List<LocalDate> patient1TraceDatesDescend = patient1TraceDates;

        for(LocalDate patient1Date : patient1TraceDatesDescend)
        {
            if(patient1Date.isAfter(onOrBefore))
                continue;

            Object patient1LocationOnDay = _patientTraces.get(patient1);
            Object patient2LocationOnDay = _patientTraces.get(patient2);

            if(patient1LocationOnDay == patient2LocationOnDay)
                return patient1Date;
        }

        throw new IllegalArgumentException("Given patients never overlap");

    }

    private LocalDate findFirstOverlap(P patient1, P patient2)  throws NoOverlapException
    {
        Map<LocalDate,Object> patient1Locations = _patientTraces.get(patient1);
        Map<LocalDate,Object> patient2Locations = _patientTraces.get(patient2);

        List<LocalDate> patient1TraceDates = new ArrayList(patient1Locations.keySet());
        Collections.sort(patient1TraceDates);
        List<LocalDate> patient1TraceDatesAscend = patient1TraceDates;

        for(LocalDate patient1Date : patient1TraceDatesAscend)
        {

            Object patient1LocationOnDay = _patientTraces.get(patient1);
            Object patient2LocationOnDay = _patientTraces.get(patient2);

            if(patient1LocationOnDay == patient2LocationOnDay)
                return patient1Date;
        }

        throw new NoOverlapException();

    }

    private SnitkinEdge<P, D> findEdge(P from, P to, Set<SnitkinEdge<P, D>> potentialResultEdges) {

        for(SnitkinEdge<P, D> edge : potentialResultEdges)
        {
            if(edge.getSource() == from && edge.getDestination() == to)
            {
                return edge;
            }
        }

        throw new IllegalArgumentException("Given edge not in set.");

    }
}

class NoOverlapException extends Exception {}
