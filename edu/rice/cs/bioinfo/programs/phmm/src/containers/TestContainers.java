package containers;

public class TestContainers {

    /**
     * @param args
     */
    public static void main(String[] args) {
        //gonna maek this array
        ListOfArrays<Integer> a = new ListOfArrays<Integer>(new IntArrayFactory(), 2);
        a.add(1);
        a.add(2);
        a.add(3);
        a.add(4);
        a.add(5);

        System.out.println(a);
        System.out.println(a.get(3));

    }

}
