package edu.rice.cs.bioinfo.library.epidemiology.transmissionMap.Snitkin2012;

import edu.rice.cs.bioinfo.library.epidemiology.transmissionMap.TransMapInferrer;
import edu.rice.cs.bioinfo.library.epidemiology.transmissionMap.TransMapResult;
import edu.rice.cs.bioinfo.library.graph.algorithms.generation.simple.CompleteDigraphFactory;
import edu.rice.cs.bioinfo.library.graph.algorithms.minimumSpanningArborescence.MinSpanArborescence;
import edu.rice.cs.bioinfo.library.graph.algorithms.minimumSpanningArborescence.MinSpanArborescenceSolver;
import edu.rice.cs.bioinfo.library.graph.algorithms.minimumSpanningArborescence.MinSpanArborescenceSolverCompleteDigraph;
import edu.rice.cs.bioinfo.library.programming.extensions.java.lang.iterable.IterableHelp;
import org.joda.time.Days;
import org.joda.time.LocalDate;
import org.joda.time.ReadablePartial;

import java.util.*;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 12/19/12
 * Time: 12:56 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class SnitkinTransMapInferrerTemplate<P,S,D extends Comparable<D>> implements TransMapInferrer<P,SnitkinEdge<P,D>>
{
    private final Set<P> _patients;

    private final Map<P, Map<LocalDate,Object>> _patientTraces;

    private final Map<P, LocalDate> _patientToFirstPositiveDate;

    private final Map<S,P> _sequencingToPatient;

    private Integer _eMax;

    public int getEMax()
    {
        if(_eMax == null)
        {
            int numDaysSpaningAllTraceDates = Days.daysBetween(IterableHelp.min(_datesInTrace), IterableHelp.max(_datesInTrace)).getDays();
            _eMax = (numDaysSpaningAllTraceDates + 1) * 2;
        }

        return _eMax;
    }

    private final Set<ReadablePartial> _datesInTrace;


     /*
    private DirectedGraph<P,SnitkinEdge<P,D>> makeEmptyGraph()
    {
        return new DirectedSparseGraph<P, SnitkinEdge<P, D>>();
    }
             */
    public SnitkinTransMapInferrerTemplate(Map<P, Map<LocalDate, Object>> patientTraces,
                                           Map<P, LocalDate> patientToFirstPositiveDate,
                                           Map<S, P> sequencingToPatient) {

        _patients = patientToFirstPositiveDate.keySet();
        _patientTraces = patientTraces;
        _patientToFirstPositiveDate = patientToFirstPositiveDate;
        _sequencingToPatient = sequencingToPatient;

        _datesInTrace = new HashSet<ReadablePartial>();
        for(P patient : patientTraces.keySet())
        {
            _datesInTrace.addAll(patientTraces.get(patient).keySet());
        }


    }

    public SnitkinTransMapInferrerResult<SnitkinEdge<P, D>> inferMaps()
    {
        return inferMaps(null);
    }

    public  SnitkinTransMapInferrerResult<SnitkinEdge<P, D>> inferMaps(P assumedRoot)
    {
        Set<SnitkinEdge<P,D>> completePatientGraphEdges = makeCompletePatientGraph(assumedRoot);

        return inferMapsFromEdges(completePatientGraphEdges);

    }

    SnitkinTransMapInferrerResult inferMapsFromEdges(Set<SnitkinEdge<P, D>> completePatientGraphEdges)
    {
        MinSpanArborescenceSolver<SnitkinEdge<P,D>,D> minSpanSolver = new MinSpanArborescenceSolverCompleteDigraph<SnitkinEdge<P,D>,D>(makeDistance(0), makeDistance(1))
        {
            @Override
            protected D addWeight(D w1, D w2) {
                return add(w1, w2);
            }

            @Override
            protected D subtractWeight(D w1, D w2) {
                return subtract(w1, w2);
            }

            @Override
            protected D getWeightOfEdge(SnitkinEdge<P, D> edge) {
                return edge.getDistance();
            }

            @Override
            protected P getSource(SnitkinEdge<P, D> edge) {
                return edge.getSource();
            }

            @Override
            protected P getDestination(SnitkinEdge<P, D> edge) {
                return edge.getDestination();
            }
        };

        MinSpanArborescence<SnitkinEdge<P,D>,D> minSpanTree = minSpanSolver.tryFindMinSpanArborescence(completePatientGraphEdges);

        int geneticDistanceAccum = 0;
        int silentColonizationAccum = 0;
        for(SnitkinEdge<P,D> edge : minSpanTree.Edges)
        {
            geneticDistanceAccum += edge.getGeneticDistance();
            silentColonizationAccum += edge.getEpidemiologicalDistance();
        }

        Set<Set<SnitkinEdge<P, D>>> minimumSpans = new HashSet<Set<SnitkinEdge<P, D>>>();
        minimumSpans.add(minSpanTree.Edges);

        return new SnitkinTransMapInferrerResult(minimumSpans, geneticDistanceAccum, silentColonizationAccum);

    }

    Set<SnitkinEdge<P,D>> makeCompletePatientGraph(final P assumedRoot)
    {
        final Map<P, Map<P, Integer>> geneticDistanceBetweenPatients = computeGeneticDistanceBetweenPatients();

        return new CompleteDigraphFactory<P,SnitkinEdge<P,D>>()
               {
                   @Override
                   public SnitkinEdge<P, D> makeEdge(P source, P destination)
                   {
                       if(assumedRoot != null && destination.equals(assumedRoot))
                       {
                          return new SnitkinEdgeBase<P, D>(source,  destination, -1, -1) {
                              public D getDistance() {
                                  return getMaxDistance();
                              }
                          };
                       }
                       else
                       {
                           int geneticDistance = geneticDistanceBetweenPatients.get(source).get(destination);
                           int epidemiologicalDistance = getEpidemiologicalDistance(source,  destination);
                           SnitkinEdge<P,D> edge = SnitkinTransMapInferrerTemplate.this.makeEdge(source, destination, geneticDistance, epidemiologicalDistance);
                           return edge;
                       }
                   }
               }.makeCompleteDigraph(_patients);
    }

    protected abstract D getMaxDistance();

    private Map<P, Map<P, Integer>> computeGeneticDistanceBetweenPatients()
    {
        Map<P,Map<P,Integer>> geneticDistanceBetweenPatients = new HashMap<P, Map<P, Integer>>();
        for(P patient : _patients)
        {
            geneticDistanceBetweenPatients.put(patient, new HashMap<P, Integer>());
        }


        for(S genome1 : _sequencingToPatient.keySet())
        {
            P patientOfGenome1 = _sequencingToPatient.get(genome1);
            Map<P,Integer> patient1Distances = geneticDistanceBetweenPatients.get(patientOfGenome1);

            for(S genome2 : _sequencingToPatient.keySet())
            {
                P patientOfGenome2 = _sequencingToPatient.get(genome2);

                if(patientOfGenome1.equals(patientOfGenome2))
                    continue;

                if(!_patients.contains(patientOfGenome1) || !_patients.contains(patientOfGenome2))
                    continue;

                int geneticDistance = getGeneticDistance(genome1, genome2);

                if(patient1Distances.containsKey(patientOfGenome2))
                {
                    int currentMin = patient1Distances.get(patientOfGenome2);
                    if(currentMin > geneticDistance)
                    {
                        patient1Distances.put(patientOfGenome2, geneticDistance);
                    }
                }
                else
                {
                    patient1Distances.put(patientOfGenome2, geneticDistance);
                }

            }
        }
        return geneticDistanceBetweenPatients;
    }

    protected int getEpidemiologicalDistance(P donor, P recipient)
    {
        LocalDate donorFirstPositive = _patientToFirstPositiveDate.get(donor); // the first day the donor tested positive
        LocalDate recipientFirstPositive = _patientToFirstPositiveDate.get(recipient);  // the first day the recipient tested positive

        // find the first day both the donor and recipient were in the same place,
        // or null if they were never in the same place.
        LocalDate firstOverlap = null;
        try
        {
            firstOverlap = findFirstOverlap(donor, recipient);
        }
        catch (NoOverlapException e)
        {
            return this.getEMax(); // no overlap.  Unlikely event, so make edge weight very large.
        }

        if(recipientFirstPositive.isBefore(firstOverlap))
        {
            return this.getEMax(); // recipient tests positive before first overlap.  Unlikely event, so make edge weight very large.
        }

        Map<LocalDate,?> donorLocations = _patientTraces.get(donor); // donor locations by day
        Map<LocalDate,?> recipientLocations = _patientTraces.get(recipient); // recipient locations by day

        LocalDate assumedTransmissionDate = null;
        for(LocalDate i = recipientFirstPositive; assumedTransmissionDate == null; i = i.minusDays(1))
        {
            boolean patientsAreCollocated;
            if(donorLocations.containsKey(i) && recipientLocations.containsKey(i)) // do we know both locations on day i?
            {
                Object donorLocation = donorLocations.get(i);
                Object recipientLocation = recipientLocations.get(i);
                patientsAreCollocated = donorLocation.equals(recipientLocation);
            }
            else
            {
                patientsAreCollocated = false;
            }

            if(patientsAreCollocated) // because firstOverlap is not null, must be true at some point in time.
            {
                assumedTransmissionDate = i;
            }
        }

        int recipSilentColonizationDaysCount = Days.daysBetween(assumedTransmissionDate, recipientFirstPositive).getDays();
        int donorSilentColonizationDaysCount = donorFirstPositive.isAfter(assumedTransmissionDate) ?
                                                Days.daysBetween(assumedTransmissionDate, donorFirstPositive).getDays() : 0;

        return recipSilentColonizationDaysCount + donorSilentColonizationDaysCount;


    }

    /*
    private D getEpidemiologicalDistance(P donor, P recipient)
    {
        LocalDate donorFirstPositive = _patientToFirstPositiveDate.get(donor);
        LocalDate recipientFirstPositive = _patientToFirstPositiveDate.get(recipient);

        LocalDate firstOverlapOnBeforeRecipFirstPositive;
        try
        {
            firstOverlapOnBeforeRecipFirstPositive =  findLastDayOfMostRecentOverlap(donor, recipient, recipientFirstPositive);
        }
        catch (NoOverlapException e)
        {
            return getEMax();
        }

        if(recipientFirstPositive.isBefore(donorFirstPositive) || !recipientFirstPositive.isAfter(donorFirstPositive))
        {
            int numDaysBetween = Days.daysBetween(firstOverlapOnBeforeRecipFirstPositive, donorFirstPositive).getDays();
            return makeDistance(numDaysBetween);

        }
        else
        {
            int numDaysBetween = Days.daysBetween(firstOverlapOnBeforeRecipFirstPositive, recipientFirstPositive).getDays();
            return makeDistance(numDaysBetween);
        }


    }  */



    private D getSilentColonizationOfRecipient(P donor, P recipient, LocalDate donorFirstPositive, LocalDate recipientFirstPositive) throws NoOverlapException {

        LocalDate firstOverlapAfterDonorPositive = findFirstDayOfMostRecentOverlap(donor, recipient, recipientFirstPositive);
        int numDays = Days.daysBetween(firstOverlapAfterDonorPositive, recipientFirstPositive).getDays();

        if(numDays < 0)
        {
            return makeDistance(0);
        }

        return makeDistance(numDays);
    }

    private D getSilentColonizationOfDonor(LocalDate donorFirstPositive, LocalDate recipientFirstPositive)
    {
        int numDays = Days.daysBetween(recipientFirstPositive, donorFirstPositive).getDays();
        if(numDays > 0)
        {
            return makeDistance(numDays);
        }
        else
        {
            return makeDistance(0);
        }
    }

    private LocalDate findFirstDayOfMostRecentOverlap(P patient1, P patient2, LocalDate onOrBefore) throws NoOverlapException {

        Map<LocalDate,?> patient1Locations = _patientTraces.get(patient1);
        Map<LocalDate,?> patient2Locations = _patientTraces.get(patient2);

        List<LocalDate> patient1TraceDatesAscend = new ArrayList<LocalDate>(patient1Locations.keySet());
        Collections.sort(patient1TraceDatesAscend);
        LocalDate dayBeforeFirstDateOfDataOnPatient1 = patient1TraceDatesAscend.get(0).minusDays(1);

        List<LocalDate> patient2TraceDatesAscend = new ArrayList<LocalDate>(patient2Locations.keySet());
        Collections.sort(patient2TraceDatesAscend);
        LocalDate dayBeforeFirstDateOfDataOnPatient2 = patient2TraceDatesAscend.get(0).minusDays(1);

        for(LocalDate i = onOrBefore; dayBeforeFirstDateOfDataOnPatient1.isBefore(i) && dayBeforeFirstDateOfDataOnPatient2.isBefore(i); i = i.minusDays(1))
        {
            if(!patient1Locations.containsKey(i) || !patient2Locations.containsKey(i))
            {
                continue;
            }

            Object patient1LocationOnI = patient1Locations.get(i);
            Object patient2LocationOnI = patient2Locations.get(i);

            if(patient1LocationOnI.equals(patient2LocationOnI))
            {
                do {
                    i = i.minusDays(1);
                    if(!patient1Locations.containsKey(i) || !patient2Locations.containsKey(i))
                    {
                        continue;
                    }
                    patient1LocationOnI = patient1Locations.get(i);
                    patient2LocationOnI = patient2Locations.get(i);
                    if(!patient1LocationOnI.equals(patient2LocationOnI))
                    {
                        return i.plusDays(1);
                    }
                }
                while (dayBeforeFirstDateOfDataOnPatient1.isBefore(i) && dayBeforeFirstDateOfDataOnPatient2.isBefore(i));

                return i.plusDays(1);
            }
        }

        throw new NoOverlapException();

    }

    private LocalDate findLastDayOfMostRecentOverlap(P patient1, P patient2, LocalDate onOrBefore) throws NoOverlapException {

            Map<LocalDate,?> patient1Locations = _patientTraces.get(patient1);
            Map<LocalDate,?> patient2Locations = _patientTraces.get(patient2);

            List<LocalDate> patient1TraceDatesAscend = new ArrayList<LocalDate>(patient1Locations.keySet());
            Collections.sort(patient1TraceDatesAscend);
            LocalDate dayBeforeFirstDateOfDataOnPatient1 = patient1TraceDatesAscend.get(0).minusDays(1);

            List<LocalDate> patient2TraceDatesAscend = new ArrayList<LocalDate>(patient2Locations.keySet());
            Collections.sort(patient2TraceDatesAscend);
            LocalDate dayBeforeFirstDateOfDataOnPatient2 = patient2TraceDatesAscend.get(0).minusDays(1);

            for(LocalDate i = onOrBefore; dayBeforeFirstDateOfDataOnPatient1.isBefore(i) && dayBeforeFirstDateOfDataOnPatient2.isBefore(i); i = i.minusDays(1))
            {
                if(!patient1Locations.containsKey(i) || !patient2Locations.containsKey(i))
                {
                    continue;
                }

                Object patient1LocationOnI = patient1Locations.get(i);
                Object patient2LocationOnI = patient2Locations.get(i);

                if(patient1LocationOnI.equals(patient2LocationOnI))
                {
                    return i;
                }
            }

            throw new NoOverlapException();

        }

    private LocalDate findFirstOverlap(P donor, P recipient)  throws NoOverlapException
    {
        return findFirstOverlap(donor, recipient, null);
    }

    private LocalDate findFirstOverlap(P donor, P recipient, LocalDate onOrAfter)  throws NoOverlapException
    {
        Map<LocalDate,?> donorLocations = _patientTraces.get(donor);
        Map<LocalDate,?> recipientLocations = _patientTraces.get(recipient);

        List<LocalDate> donorTraceDates = new ArrayList<LocalDate>(donorLocations.keySet());
        Collections.sort(donorTraceDates);
        List<LocalDate> donorTraceDatesAscend = donorTraceDates;

        for(LocalDate donorDate : donorTraceDatesAscend)
        {
            if(onOrAfter != null && donorDate.isBefore(onOrAfter))
            {
                continue;
            }

            if(recipientLocations.containsKey(donorDate))
            {
                Object donorLocationOnDay = donorLocations.get(donorDate);
                Object recipientLocationOnDay = recipientLocations.get(donorDate);

                if(donorLocationOnDay.equals(recipientLocationOnDay))
                    return donorDate;
            }
        }

        throw new NoOverlapException();

    }

    protected abstract int getGeneticDistance(S sequence1, S sequence2);

    protected abstract D add(D a, D b);

    protected abstract D subtract(D a, D b);

    protected abstract D makeDistance(int intDistance);

    protected abstract SnitkinEdge<P,D> makeEdge(P source, P destination, int geneticDistance, int epidemiologicalDistance);

}

class NoOverlapException extends Exception {}
