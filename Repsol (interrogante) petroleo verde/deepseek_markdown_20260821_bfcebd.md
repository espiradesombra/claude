# Proyecto CERO_EMISIONES_REPSOL_33x1

**Autor:** Víctor Manzanares Alberola (VMA)
**Fecha:** 21 de Agosto de 2026
**Estado:** Propuesta Técnica y Financiera para Repsol

---

## Resumen del Proyecto

Este repositorio contiene el diseño, la simulación y el caso de negocio para un sistema de **Ciclo Cerrado de Combustible** que permite a vehículos de combustión (gasoil/diésel) operar con emisiones netas casi nulas y una autonomía significativamente mejorada.

El sistema propuesto integra cuatro tecnologías clave en un solo vehículo:

1.  **Motor 3×O2:** Inyección de aire presurizado para una combustión completa y eficiente.
2.  **Captura DAC:** Extracción y almacenamiento a bordo del CO2 producido.
3.  **Síntesis Sabatier:** Reconversión del CO2 capturado en combustible sintético.
4.  **Batería Kilómetro:** Fuente de energía adicional y sostenible para alimentar el proceso de síntesis.

Este proyecto es una pieza fundamental del ecosistema **33×1**, que propone un intercambio de tecnología por 33 años de paz global.

## Estructura del Repositorio

*   **`docs/`**: Contiene la documentación técnica y regulatoria.
*   **`simulaciones/`**: Scripts en Python para validar la física y la termodinámica del sistema.
*   **`financiero/`**: Análisis de viabilidad económica y casos de uso.
*   **`presentacion/`**: Documentos para la presentación formal a Repsol.

## Primeros Pasos

Para ejecutar las simulaciones y validar los cálculos:

1.  Asegúrate de tener Python 3.8+ instalado.
2.  Instala las dependencias: `pip install -r simulaciones/requirements.txt`
3.  Ejecuta los gemelos en orden:
    ```bash
    python simulaciones/gemelo_1_dac.py
    python simulaciones/gemelo_2_sintesis.py
    python simulaciones/gemelo_3_ciclo_completo.py