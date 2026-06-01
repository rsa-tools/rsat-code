# 🛡️ PR de Seguridad – RSAT

## 🚦 Estado del PR

- [ ] Borrador / análisis inicial
- [ ] Implementación en curso
- [ ] Implementación terminada
- [ ] Probado en VM
- [ ] Pruebas funcionales completadas
- [ ] Pruebas de seguridad completadas
- [ ] Validado por Jacques
- [ ] Listo para merge


## 📌 Resumen

Describe brevemente qué problema de seguridad se corrige.

---

## 🔍 Problema identificado

- Tipo de vulnerabilidad:
  - [ ] Command Injection
  - [ ] Path Traversal
  - [ ] Arbitrary File Read
  - [ ] Input Validation
  - [ ] Otro: ___________

- Descripción:

(Explicar claramente el problema)

---

## ⚠️ Riesgo

- [ ] Ejecución remota de comandos (RCE)
- [ ] Lectura de archivos del sistema
- [ ] Exposición de datos
- [ ] Otro: ___________

---

## 🛠️ Cambios realizados

- Parámetro afectado:
- Estrategia aplicada:
- Archivos modificados:

```text
(agregar archivos)
```


## 🔐 Nueva estrategia implementada

(Explicar cómo se resuelve ahora el problema)



## 🧪 Pruebas de seguridad

### Caso: Path traversal

Input:

```text
../../../../etc/passwd
```

Resultado esperado:

```text
Acceso denegado
```

Resultado obtenido:

 ```text
(describir resultado)
 ```


### Caso: Command injection

Input:

 ```text
origin=genomic; ls
 ```

Resultado esperado:

 ```text
Rechazo del parámetro
 ```

Resultado obtenido:

 ```text
(describir resultado)
 ```

---

## ✅ Pruebas funcionales

- [ ] retrieve-seq funciona
- [ ] matrix-scan funciona
- [ ] resultados correctos

Descripción:

 ```text
(describir pruebas)
 ```

---

## 📸 Evidencia

(Logs, outputs o capturas)

---

## ⚙️ Impacto

- [ ] No rompe compatibilidad
- [ ] Cambios mínimos
- [ ] Requiere documentación

Descripción:

```text
(explicar impacto)
```

---

## 📚 Documentación relacionada

Issue:

```text
#ID
```

---

## 👥 Revisión

- Implementado por:
- Revisado por:
- Validación final:

---

## 🚀 Checklist final

- [ ] Vulnerabilidad corregida
- [ ] Inputs validados
- [ ] Sin uso de shell
- [ ] Pruebas de seguridad ejecutadas
- [ ] Pruebas funcionales ejecutadas
- [ ] Probado en VM
- [ ] Documentación actualizada
- [ ] Listo para validación

---

## 🧠 Notas adicionales

```text
(comentarios)
```
