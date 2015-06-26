package edu.rice.cs.bioinfo.programs.phylonet.algos.hmmEthan.gui;


import com.google.gson.Gson;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.data.xy.*;

import javax.swing.*;
import java.awt.*;

public class PlotBarGraph extends SavableTabData {

    JFreeChart chart;
    String title;
    String XAXISTITLE;
    String typeOfData;
    int[] inputData;

    /**
     * Deserializes a PlotAreaGraph from a JsonObject
     *
     * @param obj The JsonObject holding the data.
     */
    public PlotBarGraph(JsonObject obj) {
        this(obj.get("title").getAsString(),
                obj.get("XAXISTITLE").getAsString(),
                obj.get("typeOfData").getAsString()
                , new Gson().fromJson(obj.get("values"), int[].class));
    }

    /**
     * Creates a PlotAreaGraph from data with a given title.
     *
     * @param title      The title of the graph.
     * @param typeOfData The type of the data (IE, the y axis)
     * @param inputData  The data, where inputData[x] is the y value for a type at a given x coordinate.
     */
    public PlotBarGraph(String title, String XAXIS, String typeOfData, int[] inputData) {
        this.title = title;
        XAXISTITLE = XAXIS;
        this.typeOfData = typeOfData;
        this.inputData = inputData;

        chart = createChart(title, typeOfData, inputData);
    }


    /**
     * Creates a dataset from a two dimensional array.
     *
     * @param inputData The data, where inputData[x][type] is the y value for a type at a given x coordinate.
     * @return A dataset holding the given data.
     */
    private XYDataset createDataset(int[] inputData) {
        XYSeriesCollection collection = new XYSeriesCollection();

        collection.addSeries(createSeries("", inputData));
        return collection;
    }

    /**
     * Creates a series from data with a given label.
     *
     * @param label   The label for the series.
     * @param data The data for the series.
     * @return The corresponding XYSeries.
     */
    private XYSeries createSeries(Comparable label, int[] data) {
        XYSeries series = new XYSeries(label, false, false);
        for (int i = 0; i < data.length; i++) {
            series.add(i, data[i]);
        }
        return series;
    }

    /**
     * Creates a stacked XY area chart with the given title and dataset.
     *
     * @param title      The title of the graph.
     * @param typeOfData The type of the data (IE, the y axis)
     * @param inputdataset    The dataset to display
     * @return The chart.
     */
    private JFreeChart createChart(String title, String typeOfData, int[] inputdataset) {



        chart = ChartFactory.createXYAreaChart(
                    title, XAXISTITLE, typeOfData, createDataset(inputdataset), PlotOrientation.VERTICAL, true, true, false);


        chart.setBackgroundPaint(Color.white);

        final XYPlot plot = (XYPlot) chart.getPlot();
        plot.setBackgroundPaint(Color.lightGray);
        plot.setRangeGridlinePaint(Color.white);

        // customise the range axis...
        final NumberAxis rangeAxis = (NumberAxis) plot.getRangeAxis();
        rangeAxis.setStandardTickUnits(NumberAxis.createIntegerTickUnits());
        rangeAxis.setAutoRangeIncludesZero(true);

        final NumberAxis domainAxis = (NumberAxis) plot.getDomainAxis();
        domainAxis.setRange(0,inputdataset.length);
        return chart;
    }

    @Override
    public JComponent getComponent() {
        return new ChartPanel(chart);
    }

    @Override
    public String getType() {
        return "PlotBarGraph";
    }

    @Override
    public JsonElement getData() {

        Gson gson = new Gson();

        JsonObject obj = new JsonObject();
        obj.addProperty("title", title);
        obj.addProperty("XAXISTITLE",XAXISTITLE);
        obj.addProperty("typeOfData", typeOfData);
        obj.add("values", gson.toJsonTree(inputData));

        return obj;
    }

}