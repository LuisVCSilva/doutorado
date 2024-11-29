using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class SimManager : MonoBehaviour
{
    private float taValue = 0f;
    private float tdValue = 0f;
    public GameObject lattice;
    // Variáveis para a SelectionGrid
    private int selectedIndex = 0;
    private string[] options = { "8", "16", "32", "64", "128", "256" };

    void OnGUI()
    {
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
        GameObject.Instantiate(lattice,Vector3.zero,Quaternion.identity);
        lattice.GetComponent<Lattice>().Initialize(32,32,0.1f,1.0f);
        }

        if (GUI.Button(new Rect(100, 160, 80, 30), "Pausar"))
        {
            Debug.Log("Simulação pausada.");
        }

        if (GUI.Button(new Rect(190, 160, 80, 30), "Resetar"))
        {
            Debug.Log("Valores resetados.");
            taValue = 0;
            tdValue = 0;
            selectedIndex = 0;
        }
    }
}
