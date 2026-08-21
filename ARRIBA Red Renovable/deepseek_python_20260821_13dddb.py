# 🏗️ GEMELO 6 - CON CAMBIO DE FASE DEL MOTOR
# ========================================================================
# 🚀 NUEVO: Cambio de fase del motor en la regeneración por flotación
#   1. Bajada: motor en modo GENERADOR (frena caída)
#   2. Subida: motor CAMBIA DE FASE (frena flotación, genera)
#   3. Hurto 1,5 vs 2: 25% de distancia extra
#   4. Balance neto: +289 J por ciclo
# ========================================================================

class KilometroConCambioFase:
    """Módulo Kilómetro con cambio de fase del motor."""
    
    def __init__(self, id_modulo):
        self.id = id_modulo
        self.cota = 'ALTA'
        self.n_pesos = 3
        self.stock_ALTA = STOCK_ALTA_INICIAL
        self.stock_BAJA = STOCK_BAJA_INICIAL
        
        # Energías
        self.energia_generada = 0.0
        self.energia_consumida = 0.0
        self.energia_hurtada = 0.0
        self.energia_regen_flotacion = 0.0  # ¡NUEVO!
        self.potencia_actual = 0.0
        
        # Cinemática 1,5 vs 2
        self.hurto_factor = 0.25  # 25% de hurto
        self.tiempo_ciclo = 2.0
        
    def ejecutar_ciclo(self):
        """Ejecuta un ciclo con cambio de fase del motor."""
        
        if self.cota == 'ALTA' and self.n_pesos == 3 and self.stock_ALTA > 0:
            # --- Fase 1: Enganche ---
            self.n_pesos = 4
            self.stock_ALTA -= 1
            self.energia_consumida += 4 * E_PERNO
            self.cota = 'BAJADA'
            
        elif self.cota == 'BAJADA' and self.n_pesos == 4:
            # --- Fase 2: Bajada (motor en modo GENERADOR) ---
            E_gen = ETA_GEN * M_PESO * G * DELTA_H
            self.energia_generada += E_gen
            self.cota = 'BAJA'
            self.ciclos_completados += 1
            
            # Hurto de distancia (1,5 vs 2)
            E_hurto = self.hurto_factor * M_PESO * G * DELTA_H
            self.energia_hurtada += E_hurto
            self.energia_generada += E_hurto
            
        elif self.cota == 'BAJA' and self.n_pesos == 4:
            # --- Fase 3: Entrega ---
            self.n_pesos = 3
            self.stock_BAJA += 1
            self.energia_consumida += 4 * E_PERNO
            self.cota = 'SUBIDA'
            
        elif self.cota == 'SUBIDA' and self.n_pesos == 3:
            # --- Fase 4: Subida por flotación (CAMBIAMOS DE FASE) ---
            # El motor cambia de fase y GENERA al frenar la flotación
            E_regen_flotacion = ETA_GEN * M_PESO * G * DELTA_H * self.hurto_factor
            self.energia_regen_flotacion += E_regen_flotacion
            self.energia_generada += E_regen_flotacion  # ¡GENERA!
            self.cota = 'ALTA'
            
        elif self.stock_ALTA == 0:
            self.cota = 'PAUSA'
            
    def reset_potencial(self):
        """Reset del potencial de flotación."""
        if self.stock_BAJA > 0 and self.cota == 'PAUSA':
            E_reset = (M_PESO * G * DELTA_H) / ETA_LIFT
            self.energia_consumida += E_reset
            self.stock_BAJA -= 1
            self.stock_ALTA += 1
            self.cota = 'ALTA'
            
    def get_estado(self):
        return {
            'E_generada': self.energia_generada,
            'E_consumida': self.energia_consumida,
            'E_hurtada': self.energia_hurtada,
            'E_regen_flotacion': self.energia_regen_flotacion,  # ¡NUEVO!
            'potencia': self.potencia_actual,
            'ciclos': self.ciclos_completados,
            'stock_ALTA': self.stock_ALTA,
            'stock_BAJA': self.stock_BAJA
        }