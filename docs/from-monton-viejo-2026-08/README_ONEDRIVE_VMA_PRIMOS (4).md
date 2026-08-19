# Carpeta OneDrive — Proyecto Primos (VMA)

**Fecha:** 2025-09-09

## Estructura sugerida

```
/Primos_VMA/
  ├─ Documento/
  │   ├─ VMA_Primos_Completo_v3.docx
  │   └─ PDF/ (exportados desde Word)
  ├─ Figuras/
  │   ├─ fig1_region_triangulo.png
  │   ├─ fig2_Lmenosm.png
  │   ├─ fig3_criba_hibrida.png
  │   └─ fig4_restos_scatter.png
  ├─ Datos/
  │   └─ Lm_dataset_v3.csv
  └─ Scripts/
      └─ anexoE_L_m_script.py
```

## Pasos recomendados
1. **Sube** estos archivos a la carpeta OneDrive con la estructura sugerida.
2. Abre `VMA_Primos_Completo_v3.docx` y **actualiza la tabla de contenidos**: *Referencias → Tabla de contenidos → Actualizar toda la tabla*.
3. Exporta PDF desde Word y guárdalo en `Documento/PDF/`.
4. Crea un **enlace de uso compartido** en OneDrive (Lectura) del `.docx` y del `.pdf`.
5. En tu blog (WordPress), inserta el botón/enlace con el siguiente snippet HTML (sustituye las URL por las reales de OneDrive):

```html
<p><a class="wp-block-button__link" href="https://1drv.ms/u/s!URL_DOCX" target="_blank" rel="noopener">📄 Descargar DOCX — Cribas y cotas (v3)</a></p>
<p><a class="wp-block-button__link" href="https://1drv.ms/u/s!URL_PDF" target="_blank" rel="noopener">🔗 Versión PDF — para lectura</a></p>
```

## Notas
- `Lm_dataset_v3.csv` contiene `n, L(n), K, m_sum(n), L_minus_m` para un rango amplio de n.
- `anexoE_L_m_script.py` reproduce los cálculos de L, m y L−m mostrados en el documento.
- Las figuras PNG están optimizadas para Word/Blog.
