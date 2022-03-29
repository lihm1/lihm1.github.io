---
title: hello hexo
date: 2022-03-29 15:23:35
index_img: /img/download3.jpg 
categories:
- blog
tags:
- hexo
---

参考
https://fluid-dev.github.io/hexo-fluid-docs/guide/#%E5%88%86%E7%B1%BB%E9%A1%B5
https://blog.csdn.net/yaorongke/article/details/119089190
### 1.安装Node.js
官网下载最新版
https://nodejs.org/en/
最后安装好之后,打开命令提示符，输入node -v和npm -v，出现版本号就成功了。
### 2.添加国内镜像源
如果没有梯子的话，可以使用阿里的国内镜像进行加速。
npm config set registry https://registry.npm.taobao.org
### 3.安装hexo
输入npm install -g hexo-cli  安装hexo
输入 hexo -h 检查是否安装成功
### 4.新建文件夹blog
进入我们新建的文件夹，初始化文件夹 hexo init
{% asset_img one.png This is an example image %}
我们会得到相应的hexo搭建博客的文件：
{% asset_img two.png This is an example image %}
### 5.新建博客
hexo new post "your blog name"
打开文件夹的source/_posts / 里面包含自己的博客的md文件
### 6.生成博客静态网址
~~~
hexo g 
hexo s
~~~
{% asset_img three.png This is an example image %}
{% asset_img four.png This is an example image %}
输入网址就可以看到自己的博客了
### 7.博客的主题修改
搜索自己喜欢的主题，将主题下载到blog/themes 文件夹下。并将文件名修改为主题名。
下面以fluid的主题为例：
最新版本下载地址为：https://github.com/fluid-dev/hexo-theme-fluid/releases
下载后解压，并重命名为fluid。
在hexo/blog/_config.yml 修改自己的主题
#Extensions
##Plugins: https://hexo.io/plugins/
##Themes: https://hexo.io/themes/
theme: fluid

查看博客能否生成：
~~~
hexo s
~~~
此时博客的样子：
{% asset_img five.png This is an example image %}
### 8.个性化博客主题
#### 8.1写文章的时候生成文件夹保存图片
将_config.yml里面的post_asset_folder改为true。
post_asset_folder: true
~~~
hexo new post "test"
~~~
{% asset_img six.png This is an example image %}
我们可以下载一张图片放入test文件夹里面，如何在文中引用，请参照https://hexo.io/zh-cn/docs/asset-folders.html
在test.md中加入下载的图片其中一种方法如下所示：
~~~
---
title: test
date: 2022-03-29 15:23:35
tags:
---
测试图片
{% asset_img download.jpg This is an example image %}
~~~
提交修改后再次打开我们的博客可以看到：

{% asset_img seven.jpg This is an example image %}
 #### 8.2文章在首页的封面图
{% asset_img eight.png This is an example image %}
设置方式如下：在文章中加上index_img
并且将图片保存在/themes/fluid/source/img中,在文章中只需要指定在img中就可以。也可以使用外链 Url 的绝对路径。
---
title: test
date: 2022-03-29 15:23:35
index_img: /img/download.jpg 
tags:
---
 #### 8.3 文章顶部放大图：
默认显示主题配置中的 post.banner_img，如需要设置单个文章的 Banner，在 文档中指定 banner_img 属性。图片仍然放置在/themes/fluid/source/img
---
title: 文章标题
tags: [Hexo, Fluid]
index_img: /img/example.jpg
banner_img: /img/post_banner.jpg
date: 2019-10-10 10:00:00
---

 #### 8.4 添加tag  categories
title: test
date: 2022-03-29 15:23:35
index_img: /img/download.jpg 
banner_img: /img/download.jpg 
categories:
- bioinformation
tags:
- software

### 9.博客与git hub 链接，并且发布到GitHub上
创建自己的GitHub账号，登陆后新建仓库命名为lihm1.github.io
在仓库中新建文件index.html，里面随便写点内容，
然后打开上方的setting,找到pages，我们可以看到自己的网址
{% asset_img nine.png This is an example image %}
发布到GitHub添加自己的github的地址以及自己的博客存在什么分支
安装hexo-deployer-git
npm install hexo-deployer-git --save
修改_config.yml的内容为自己的GitHub地址
deploy:
  type: git
  repository: https://github.com/lihm1/lihm1.github.io.git
  branch: master

将内容更新后上传到github
hexo g
hexo d
会弹出页面让登录github并且输入密码。我们再次打开自己的github就可以看到自己的仓库已经被更新。

### 10.博客源代码的保留
github创建一个分支，将自己的源码保存在分支中，以后换电脑想要更新博客的时候可以使用。
我们可以先将github上的博客的仓库克隆下来
git clone git@github.com:lihm1/lihm1.github.io.git
进入本地的lihm1.github.io.git文件夹，创建分支source
git branch source
git branch
切换到分支source
git checkout source
我们可以在本地先把lihm1.github.io.git 文件全部删掉，然后将开始新建的blog文件夹的内容复制过来。后面将修改内容以及创建的分支push 到远程仓库。
git add .
git commit -m "hexo"
git push orgin source:source
后面如果换电脑，我们可以先从github上面clone下来source分支仓库，最好设置source为默认的仓库，在setting里面可以设置。这样code的代码就是source的。
下载下来后我们需要，然后安装、初始化hexo，重新设置hexo环境。
~~~
npm install hexo
npm install hexo-deployer-git -save
~~~
然后hexo s,便可以看到自己的博客了，如果更新了内容并且想提交到github上，与原来是一样的操作参考9.发布到GitHub的内容。



