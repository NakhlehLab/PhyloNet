package edu.rice.cs.bioinfo.library.programming;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 10/24/12
 * Time: 3:13 PM
 * To change this template use File | Settings | File Templates.
 */
public class Func1ToString<T> implements Func1<T,String>
{

    public String execute(T input) {
        return input.toString();
    }
}
