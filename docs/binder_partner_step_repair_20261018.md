# 六复合物受体页面范围修复（2026-10-18）

用户指定的v2解压站点复现：AIB234第四步打开受体后切到第五步，state.step=5、showPartner=true，Protein partner组件仍可见。原开关验收只覆盖第四步内部，没有覆盖跨页面的显示约束。

web/app.js的setStep现在离开第四步立即清除showPartner并隐藏受体，同步按钮。开关回调仅在第四步执行，异步失败也不能在其他步骤恢复showPartner。沿用同一viewer和独立链组件，不修改结构、突变选择或科学评分。

验收范围：六个复合物4→1/2/3/5/6/7、返回4、突变定位、加入/清空空间比较、cartoon/surface切换、重置相机、同一viewer、file离线零HTTP与零pageerror。14项生命周期/离线导出回归通过（8.20秒）；JS语法及按cr-at-eol检查diff通过。

证据根目录workspace_data/20261016_raw data_binder/partner_step_repair_20261018。browser_acceptance.json逐例记录，frontend_changes.json逐文件绑定修复前后SHA，delivery.json只有六例浏览器通过且最终ZIP CRC/成员/前端字节检查通过才写PASS。原v2 ZIP及修复前脚本保留；实际更新用户指定解压目录，另存v3分享ZIP。

未重算模型、未跑全仓库测试、未重新实测在线页面或五个未改AF3页面。完整交付以delivery.json为准。

## 定位结构后表示重建复发（2026-10-18）

用户AIB193（84aa、S1D/R4A）截图证明此前修复与验收不完整。真实浏览器记录：定位前整体Polymer及Protein partner父/子表示均隐藏；定位后父组件仍hidden=true，但新representation.hidden=false、repr.visible=true、renderObjects均true。根因是PDBe Mol* applyVisualStyles重建表示未继承隐藏状态，旧父组件检查产生假阳性。

renderer.js现用applyChainVisibility在表示重建、着色及开关后统一隐藏原完整结构子树、同步partner子树；重载回退重建链组件。app.js删除selectDesign在setRepresentation之前重复启动的repaintViews，避免并发着色与表示删除造成state ref错误。

15项生命周期/离线导出回归通过（4.52秒），新表示可见性和先重建再着色两项检查均使旧代码失败。验收脚本现核对实际表示/渲染对象可见性和可见原子链集合，捕获viewerError；另保存六例定位后截图。证据根目录为workspace_data/20261016_raw data_binder/partner_representation_repair_20261018，最终browser_acceptance.json及delivery.json为本次终态；上次partner_step_repair的PASS只覆盖父组件状态，不代表解决此次缺陷。

六例仅更新前端assets及其manifest哈希，完整逐文件before/after清单保留。原v2/v3 ZIP不改，新v4分享包包含修复。最终交付以delivery.json=PASS为准；未重算模型、未重验在线版或未改的五例AF3。

六例新浏览器验收全部通过：AIB193/AIB234/AIB226/AIB239/AIB201/AIB15定位后实际可见组件均仅Binder、可见链均仅A；其他六步与界面开关、多选/清空、表示切换保持；viewerError/pageerror/HTTP均零。逐例完整可见性快照与截图已写入本次browser_acceptance.json及同目录PNG。归档终态继续以delivery.json为准。
