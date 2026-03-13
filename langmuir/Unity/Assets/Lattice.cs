using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Lattice : MonoBehaviour
{
    public GameObject cubePrefab;
    public int Lx = 0;
    public int Ly = 0;

    public double rA = 0.0;
    public double rD = 0.0;
    public float spacing = 1.0f;
    private GameObject[,] grid;

    private int totalIterations;
    public float condD, condA;

    private float t = 0f;
    private int ic = 0;
    public int contagem = 0;
    private List<float> historicoTheta;
    private List<float> historicoTempo;
    private float area;
    public float rc = 0.0f;

    private float taValue = 0f;
    private float tdValue = 0f;
    // Variáveis para a SelectionGrid
    private int selectedIndex = 0;
    private string[] options = { "8", "16", "32", "64", "128", "256" };
    public bool isSimulating = false;



    void OnGUI() {
        // Definir estilo básico para os rótulos
        GUIStyle labelStyle = new GUIStyle(GUI.skin.label);
        labelStyle.fontSize = 14;

        // Slider "Ta" com rótulo
        GUI.Label(new Rect(10, 10, 50, 20), "Ta", labelStyle);
        taValue = GUI.HorizontalSlider(new Rect(60, 15, Screen.width/4, 20), taValue, 0, 10);
        GUI.Label(new Rect(270, 10, 50, 20), taValue.ToString("F2"), labelStyle);

        // Slider "Td" com rótulo
        GUI.Label(new Rect(10, 50, 50, 20), "Td", labelStyle);
        tdValue = GUI.HorizontalSlider(new Rect(60, 55, Screen.width/4, 20), tdValue, 0, 10);
        GUI.Label(new Rect(270, 50, 50, 20), tdValue.ToString("F2"), labelStyle);

        // SelectionGrid para seleção de números
        GUI.Label(new Rect(10, 90, 200, 20), "Res. Lattice:", labelStyle);
        selectedIndex = GUI.SelectionGrid(new Rect(10, 110, 300, 30), selectedIndex, options, options.Length);

        // Botões "Simular", "Pausar" e "Resetar"
        if (GUI.Button(new Rect(10, 160, 80, 30), "Simular")) {
            int res = selectedIndex==0 ? 8 : (selectedIndex==1 ? 16 : selectedIndex==2 ? 32 : selectedIndex==3 ? 64 : selectedIndex==4 ? 64 : selectedIndex==5 ? 128 : 256);
            Lx = res;
            Ly = res;
            rA = taValue;
            rD = tdValue;
            Vector3 centroid = Initialize(res,res,taValue,tdValue);
            isSimulating = true;
            Camera.main.transform.position = new Vector3(centroid.x,centroid.y,-10);
            Camera.main.orthographicSize = Lx;
        }

        if (GUI.Button(new Rect(100, 160, 80, 30), isSimulating ? "Pausar" : "Play"))
        {
            isSimulating = !isSimulating;
            Debug.Log("Simulação pausada.");
        }

        if (GUI.Button(new Rect(190, 160, 80, 30), "Resetar"))
        {
            isSimulating = false;
            foreach (Transform child in transform) {
            Destroy(child.gameObject);  // Destroy each child gameObject
            }
            Debug.Log("Valores resetados.");
            taValue = 0;
            tdValue = 0;
            selectedIndex = 0;
            int res = -1;
            Lx = 0;
            Ly = 0;
            rA = 0.0f;
            rD = 0.0f;
            rc = 0.0f;
            ic = 0;
            area = 0;
            t = 0;
        }

        GUI.Label(new Rect(10, Screen.height-75, 200, 20), "Tempo (s): " + t.ToString("F2"), labelStyle);
        GUI.Label(new Rect(10, Screen.height-50, 200, 20), "Cobertura (rc) aprox: " + (rc * 100f).ToString("F2") + "%", labelStyle);
        GUI.Label(new Rect(10, Screen.height-25, 200, 20), "Cobertura (rc) exata: " + (taValue/(taValue+tdValue) * 100f).ToString("F2") + "%", labelStyle);

    }

    // Initialization method
    public Vector3 Initialize(int Lx, int Ly, float rA, float rD)
    {
        Vector3 centroid = Vector3.zero;
        area = Lx * Ly;
        historicoTheta = new List<float>();
        historicoTempo = new List<float>();
        int LT = 100;



        totalIterations = LT;  // Just for the sake of simplicity. You can modify this logic as needed.

        Camera.main.orthographic = true;
        grid = new GameObject[Lx, Ly];

        for (int i = 0; i < Lx; i++)
        {
            for (int j = 0; j < Ly; j++)
            {
                condD = (float)(tdValue / (taValue + tdValue));
                condA = (float)(taValue / (taValue + tdValue));
                float r = Random.value;
                Vector3 position = transform.position + new Vector3(i * spacing, j * spacing, 0.0f);
                grid[i, j] = Instantiate(cubePrefab, position, Quaternion.identity);
                grid[i, j].transform.parent = transform;
                if(i==j && i==Lx/2) {
                    centroid = grid[i,j].transform.position;
                }
                Square square = grid[i, j].GetComponent<Square>();
                if (r < rc)
                {
                    square.Occupy();
                    ic++;
                }
                else
                {
                    square.Clear();
                }
            }
        }   
        return centroid;
    }

    // Update method for simulation
    private void Update()
    {
        if(isSimulating==true) {
        // Only perform iterations if the condition is met
        for (int _ = 0; _ < totalIterations; _++)
        {
            int i = Random.Range(0, Lx);
            int j = Random.Range(0, Ly);
            float r = Random.value;
            Square square = grid[i, j].GetComponent<Square>();
            float dt = 0.0f;
            
            // If square is occupied, handle transition
            if (square.state)
            {
                if (r <= condD)
                {
                    square.state = false;  // Clear square
                    ic--;
                    r = Random.value;  // Recalculate random value for dt
                    dt = (float)(-Mathf.Log(r) / (ic * rD + (area - ic) * rA));
                    t += dt;
                }
            }
            // If square is not occupied, handle transition
            else
            {
                if (r <= condA)
                {
                    square.state = true;  // Occupy square
                    ic++;
                    r = Random.value;  // Recalculate random value for dt
                    dt = (float)(-Mathf.Log(r) / (ic * rD + (area - ic) * rA));
                    t += dt;
                }
            }

            // Update the density ratio
            rc = ic / area;
            contagem++;
        }
    }
    }
}
