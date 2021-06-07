package edu.rice.cs.bioinfo.programs.soranus.viewModels;

import edu.rice.cs.bioinfo.library.programming.Proc1;

import java.util.HashSet;
import java.util.Set;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 1/30/13
 * Time: 6:11 PM
 * To change this template use File | Settings | File Templates.
 */
public class WorkspaceVM<DR extends DataRecord,
                         AR extends AnalysisRecord>
{

    public final String Title = "Soranus";

    private Set<DR> _sequencingsDataRecords = new HashSet<DR>();

    private Set<DR> _traceDataRecords = new HashSet<DR>();

    private Set<DR> _firstPositiveDataRecords = new HashSet<DR>();

    private Set<DR> _vaalOutDataRecords = new HashSet<DR>();

    private Set<Proc1<DR>> _dataRecordAddedListeners = new HashSet<Proc1<DR>>();

    public void addDataRecordAddedListener(Proc1<DR> listener)
    {
        _dataRecordAddedListeners.add(listener);
    }

    private Set<Proc1<AR>> _anaylisisAddedListeners = new HashSet<Proc1<AR>>();

    public void addAnalysisAddedListener(Proc1<AR> listener)
    {
        _anaylisisAddedListeners.add(listener);
    }

    private DocumentVM _focusDocument = null;

    public void setFocusDocument(DocumentVM focusDocument)
    {
        _focusDocument = focusDocument;

        for(Proc1<DocumentVM> listener : _focusDocumentChangedListeners)
        {
            listener.execute(_focusDocument);
        }
    }


    private Set<Proc1<DocumentVM>> _focusDocumentChangedListeners = new HashSet<Proc1<DocumentVM>>();

    public void addFocusDocumentChangedListener(Proc1<DocumentVM> listener)
    {
        _focusDocumentChangedListeners.add(listener);
    }


    public void addSequencingsData(DR data)
    {
         addDataHelp(data, _sequencingsDataRecords);
    }

    public void addTraceData(DR data)
    {
         addDataHelp(data, _traceDataRecords);
    }

    public void addFirstPositiveData(DR data) {
         addDataHelp(data, _firstPositiveDataRecords);
    }

    public void addVAALOutData(DR data)
    {
         addDataHelp(data, _vaalOutDataRecords);
    }


    private void addDataHelp(DR data, Set<DR> toAdd)
    {
        toAdd.add(data);
        for(Proc1<DR> dataRecordAddedListener : _dataRecordAddedListeners)
        {
            dataRecordAddedListener.execute(data);
        }
    }

    public void addAnalysis(AR analysisRecord)
    {
        for(Proc1<AR> listener : _anaylisisAddedListeners)
        {
            listener.execute(analysisRecord);
        }

    }

    public boolean isDataRecordTitleRepresentingSequencingsData(DR dataRecord)
    {
       return _sequencingsDataRecords.contains(dataRecord);
    }

    public boolean isDataRecordTitleRepresentingTraceData(DR dataRecord)
    {
        return _traceDataRecords.contains(dataRecord);
    }

    public boolean isDataRecordTitleRepresentingFirstPositiveData(DR dataRecord)
    {
       return _firstPositiveDataRecords.contains(dataRecord);
    }

}
