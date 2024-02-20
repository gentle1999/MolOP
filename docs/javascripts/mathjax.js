/*
 * @Author: TMJ
 * @Date: 2024-02-20 09:14:10
 * @LastEditors: TMJ
 * @LastEditTime: 2024-02-20 09:14:12
 * @Description: 请填写简介
 */
window.MathJax = {
    tex: {
      inlineMath: [["\\(", "\\)"]],
      displayMath: [["\\[", "\\]"]],
      processEscapes: true,
      processEnvironments: true
    },
    options: {
      ignoreHtmlClass: ".*|",
      processHtmlClass: "arithmatex"
    }
  };
  
  document$.subscribe(() => { 
    MathJax.startup.output.clearCache()
    MathJax.typesetClear()
    MathJax.texReset()
    MathJax.typesetPromise()
  })