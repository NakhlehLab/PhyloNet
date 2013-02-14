package edu.rice.cs.bioinfo.programs.soranus.views;

import edu.rice.cs.bioinfo.library.programming.Proc1;
import edu.rice.cs.bioinfo.library.programming.Proc3;
import edu.rice.cs.bioinfo.programs.soranus.viewModels.WorkspaceVM;

import java.io.File;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 2/6/13
 * Time: 4:02 PM
 * To change this template use File | Settings | File Templates.
 */
public interface WorkspaceView
{
    public void addDataAddRequestListener(Proc1<File> listener);

    public void addSnitkinTransMapAnalysisRequestedListener(Proc3<String,String,String> listener);

    public void addDataRecordSelectedListener(Proc1<WorkspaceVM.DataRecord> listener);

    public void addNeighborJoiningAnalysisRequestedListener(Proc1<String> listener);

    public void addAnalysisRecordSelectedListener(Proc1<WorkspaceVM.AnalysisRecord> listener);

    public void startView();

}
