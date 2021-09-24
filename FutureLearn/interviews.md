# Interview with Vivek Iyer

blob:https://www.futurelearn.com/27f612f4-d21b-4a9d-bc87-600166a67d53

Hello, and welcome to week two. Today, we’re going to be talking to Vivek Iyer, the head of human genetics informatics at the Wellcome Sanger Institute. Hi, Vivek. Hi there. So let’s start with the nice one. So how would you define bioinformatics? It’s very broad field, and it means many different things to different people. I think that it spans a range of capacities, OK? So really, what’s happened is this. The analysis of genetic data has become the province of high-performance compute. And not just high-performance compute, but the downstream analysis is actually, itself, a sort of a specialised discipline and itself has a bunch of specialised toolkits.

And so translating, taking and translating a biological experiment based on sequencing, for instance, now requires a chain of quite specialised individuals, all the way from when the DNA actually leaves or is prepared. And those specialised individuals span of this long pipeline from when the DNA is actually sequenced, all the way to when a scientist will actually evaluate plots and think about papers and things. And a bioinformatician will handle a lot of the computational tasks at any one of those parts of the pipeline. And you can see why it’s– when you look at it that way, it’s a very wide range of skill sets that, typically, one person doesn’t have.
They’ll actually have a skill set that will apply to one particular part of that pipeline. They may be dealing with the software that processes the flow cell results as they come off the sequencer, or they may be dealing with– all the way on the other end, they may be dealing with applying statistical models that will best tease out a signal between two populations just before you end up with– they’ll be talking to the scientist in the labs with plots and things, and everything in-between– control of high-performance clusters. So there’s no one good answer to that, because the process is quite spread out.

So with that variety of components or fields within bioinformatics itself, where did you first start out with bioinformatics? Well, I sort of fell into it. So I’m a physicist by training. My background is actually in black hole physics, and I got a PhD in that from the University of Chicago in a very specific niche. But I didn’t really want to continue as a theoretical physicist after I’d finished the degree, so I started working as a software consultant in a software consultancy in Chicago. And when we came to the UK almost 20 years ago, the market was dead for commercial Java developers.

But the Sanger Institute, there was this place called the Sanger and this person called Ewan Birney who needed a Java developer to work on an annotation tool called Apollo, which is a genomic annotation tool which is still in strong use– so it’s changed beyond recognition. So I started as a Java developer, and I had to get to grips with the Ensembl API, and I had to get to grips with just genomics, genome structure in general. So that’s how I started– it was writing Java tools to allow– and mosquito and fly annotators to annotate fly and mosquito genomes.

And then, a short time after that– maybe a year after that– it occurred to me that I wanted to do what all the people in the office around me were doing, which was building Ensembl gene models. And so I put down the Java, to a large part, and picked up and started coding in Perl, and started to deal with high-performance compute and informatic pipelines. And at the same time, I had a very strong grounding in gene structure and gene annotation from working around a whole bunch of Ensembl developers all the time. So that’s how I got into it. So moving forward from Ensembl, what do you see the role of bioinformatics playing in human genetics today?

So I run a group. We are a core informatics– we are, what do you call it, a team that supports the entire programme. So we see a lot of different aspects of the programme’s work. And, of course, we rest atop a huge amount of work in gene structure and annotation that’s actually been done before us. There are two main areas you can– we actually come in and where informatics actually does its thing, I guess. One of them is the running of almost-productionized pipelines at a very high volume. And that’s not a controversial area. And it’s certainly been the way that– so the classic is variant calling, for instance, right?

And the reason that informaticians have earned their living for a very long time is because the sets get progressively larger– a few thousand samples or exomes was large at one point. Then, a few thousand genomes became big. And now, tens of thousands of genomes is sort of moderate scale. So these things have been accelerating, but at any given point, the amount of the volume of the data that needs to be processed is typically just at the edge of what’s comfortable for somebody with– it forces people to actually have specialised teams that actually know how to drive compute to actually allow them to process all that data in good time.

And that’s been the classic situation, and it’s been that way for years, OK? So that’s not a controversial thing. What’s changed, I guess, in that is that large cohorts have been processed by and have been– national-level cohorts have actually been driven and provided almost as free resources. So the individual programmes don’t have– they don’t have a need to actually create cohorts themselves. So the emphasis has been shifting now from the process of big cohorts to the ingestion and the analysis of big cohorts alongside the smaller cohorts.

And that’s a whole different ballgame, because the tools change, the computer architectures, potentially, change, there are a number of very large cloud providers in the mix, and third party commercial providers in the mix, that leverage those cloud services. So a whole bunch of things are changing away from that classical model. And so allowing people to cross-analyse large cohorts is possibly a big shifting, a big part of the landscape and where we going to be. That’s one part. And I guess the other part is– yeah, I mean, there’s the variants, the human variation, and there’s the functional aspects of it.

And I think that in human genetics, one strong shift in the last couple of years has been a shift away from the analysis of pure human variation towards the cross-analysis of human variation and expression, for instance, at the same time– basically, functional outputs like expression, right? So our skill set has been changing to adapt to the change in the interests of the scientists at the same time. So given those changes in scientific interest, what would you say your typical day is? So some of the highlights, some of the challenges that you might face? Well, let’s see. So my typical day is actually just an incessant stream of emails, presentations, and Slack pinging back and forth.

So I’m not a particularly typical thing. But I can tell you kind of what used to happen, and certainly around my team who are actually doing the work. So right, some of them are developing pipelines, and they’re developing it in a variety of different toolkits for different purposes. So case in point, one of them is actually borrowing and adapting a single-cell RNA pipeline from a different lab, OK? So we are overseeing a sequence of labs- “overseeing” is the wrong term. We help a set of labs. Sometimes, as labs work with disparate– they have disparate people that don’t talk to each other as well as they might, we noticed that this lab has got technology that this lab can actually use.

So we are taking the technology from the side, reading as best we can and talking with all everybody to work out what these people need, and adopting the technology so that it comes over to these people. And so in that case, they’re adopting a single-cell RNA pipeline, and the toolkit that they’re using is a workflow engine called NextFlow. And the workflow itself drives a number of standard single-cell RNA tools like scrublet that take out the multiplets, and then that does PCA, and– standard stuff, right? Except that it’s being put together and stitched together as part of a NextFlow pipeline, and we’re adapting that pipeline from one lab to another. That’s part one.

The second kind of tool kit might be we have variant annotation databases that we’re exploring as an alternative to file-based human genome variation manipulations. Instead of sort of taking VCF, to VCF, to VCF, we’re exploring a database where you pile lots of human variation into the same database and then analyse it after the event. That’s an experimental approach. We have other kinds of pipelines and different toolkits for human variation that we do. And then, we have a strong administrative layer and an administrative approach that we have to prosecute, I guess is the term. So the reason I’m saying that is the following.

Human genetics sits atop, let’s just say, 4 to 5 petabytes of data, maybe 10% to 20% of which is in active use, and we have to police all of that. Because if you don’t, people tend to misuse, inadvertently misuse the data, the discs that they have– the high-performance discs that they have. So we’ve got– we’re ploughing a lot of development work into helping people keep track of their data, archive it, and sweep away the data they don’t really need. So we’ve got the information end of things with the pipelines, and the development dealing with the administrative end of things as well.

And so all of that’s happening at any given day, and any one developer will typically be receiving activity in one or the other area. Have you found that COVID helps or hinders with some of that? Has there been some benefits to some of the more remote working that we’ve been doing, and some things which might need a little more thought or change?

From where I stand, I think the net effect is negative.

I think positives of this, right, developers can focus. And as a team, we lean towards wishing to be left alone to do our jobs. So that’s actually a good thing. I think that the potential for mistakes increase– for misdirection, for miscommunication increase. And now more than ever, our role as a team within the human genetics department needs better communication and more frequent feedback than we are really able to effectively pursue by being spread out like this. So I would say the net effect is negative, but not crushingly so. I know where I’d rather be, but I imagine that some members of my team would probably disagree with me.

So given the current climate and the current circumstances, to someone just starting out in bioinformatics, or even yourself when you first started out, what would be the one piece of advice that you would give to them? Well, OK, so take your time and learn things thoroughly. So what will happen is that there will be this constant balance between having to get stuff done and learning things. And resist the temptation to learn as much as possible, to rush through and only half-understand what it is you’re supposed to understand before actually applying a toolkit to what you need to understand. So yeah, I mean, so that might mean–

and I can think of a variety of different contexts in which that could apply. But I mean, I don’t know. It can be learning a toolkit, or learning how to use a tool well. And that tool might be anything from a language to even just a specification, right? So I mean, read the VCF spec. I mean, really read the VCF spec, because if you don’t, if you just skitter your way through this, you’re always going to be sort of second-guessing what’s going to– what you think is. And by the time you actually need to know the answer, you will have lost the opportunity to go back and look, right?

And by the way, it will be too embarrassing to check with somebody in an office or something. So just that’s what I would say. Just take the time to be a bit more thorough at the beginning. It will actually pay off later. That’s what I would say. I genuinely think that’s one of the best pieces of advice I’ve heard and I wish I’d heard when I first started as well. Yeah, well, same, right. Possibly, the second thing is apply the right tool to the problem. So case in point, I think I was attempting to manipulate rectangular data with Perl for a solid year, double-looping through matrices, before I realised that R was out there.

And then, after I realised that R was out there, it took me– again, learn something thoroughly– it took me, possibly, six months of driving R badly before I worked out how to drive it well, OK? And then, that was something that was my mistake. So Stack Overflow doesn’t help with everything. So just be aware that there are better there might be better tools for the thing that you’re trying to do than what you’re actually trying to do at the time, and that was one possible example. Yeah, I think as well, definitely, asking for help. So we’re all friendly folk, and asking for help and asking questions you cannot answer by yourself is absolutely crucial.

No, in a way– in a, way, actually you were saying, are there any advantages or disadvantages to the COVID situation? I must admit, it’s pushed us into Slack in a big way, into these shared messaging platforms. And that is actually– inasmuch as it’s democratised the communication, I think it’s improved the communication between members of the department in some really weird ways, right? Because actually, when you ask for help, you don’t turn around. You can ask your friend, turn and ask your friend. But at the same time, you can actually post a message into a common channel. And then, generally, crowdsourcing for help that way is efficient for many different reasons.

You get an answer quickly, but also, people realise, oh, they’re doing that, OK, and may chip in advice which will say, well, OK, I’m glad you’re doing that, or you’ve thought about it, but maybe you want to try this, right? So in that sense, the democratisation through these channels has really been accelerated by COVID, which is a good thing. yeah. And I think that brings us to a good place to wrap this up. Thank you so much for joining us today, Vivek, and we really hope that if you have any comments or questions, you’ll leave them in the comments section below. And thank you for listening.

## R
blob:https://www.futurelearn.com/17e2b77e-87cb-40f2-ad78-66d497cf207d

[00:00:08.59] Hi. My name is Fatma Guerfali, and I am a researcher at Institut Pasteur de Tunis
in Tunisia. And I'm also a trainer, in Tunisia and abroad, for NGS data analysis and visualisation.
[00:00:23.23] I would like to discuss today the importance of data visualisation, especially in our
modern societies where we know that about 90% of the information that is transmitted to the
brain is visual. And today in fact, visual representations are so important that many studies
showed that, on average, a person would remember 80% of what she sees but only 20% of what
she reads. So we have the immense privilege to discuss today this particular topic on the
importance of data visualisation with Andy Kirk.
[00:01:03.75] Andy, you are a freelance data visualisation specialist. You're based in Yorkshire
UK. You work as a data visualisation consultant, a training provider, a teacher, an author,
speaker, and researcher. But most importantly, you're the editor of the award-winning amazing
website that is called visualizingdata.com. So in other words, this means you have an amazing
expertise in data visualisation. So Andy, can you tell us a bit more about yourself? And most of
all when did you sense that it was important to dedicate an entire career to explain and teach the
importance and power of data visualisation?
[00:01:50.59] Thank you, Fatma. Yes. So when I was a kid at school, I loved art and maths. And
I was always looking for a career pathway that would take me into a realm or a job or a vocation
that would bring these two worlds together. Now initially, I thought that maybe architecture was
the route to do that. However, that didn't kind of work out. My physics was a little bit too weak.
When I finished University, I'd worked on a number of modules around statistics and data
analysis. And so I took my career down that direction-- analysis, performance analysis, business
analysis.
[00:02:29.89] But it wasn't until February 2007 when I found this subject on the web-- data
visualisation. I didn't have any idea it was a thing and I found it. But when I found it, I realised
that the charts and the presentations that I was producing at work were insufficient, were not
good enough. But also I just had this immediate passion a lot because I recognised there was a
subject, there was discourse, there was evidence about good and bad practise.
[00:03:01.78] And so at the time, I was working at a University, and I had the opportunity to do a
master's degree through self-directed research. So long story short, I did a master's degree around
data visualisation, completely reinforced my passion for it. And then I thought, well, maybe I
should think about doing this for a career. So I went freelance, part time initially. And thankfully
over the last decade, the field has grown. And I've been very fortunate to be in position to take
advantage of that and grow with it.
[00:03:36.96] That's very nice, Andy. Thanks for sharing that with us. So I now understand
better what your website is presenting so much examples and resources on data visualisation.
And your website rapidly grew to become a really popular reference for followers in the field. So 
based on your experience in that particular field, what would be the basic principles for an
efficient data visualisation for you? And how can we envision data visualisation and use it as one
of the most compelling, persuasive, and potent practises to present even complex, large-scale
data today?
[00:04:17.91] That's right. And I think it's important that you mention the word complex there
because the universal goal of effective visualisation is to clarify, not necessarily simplify. There
are many topics, especially within science, where the fundamental topic is complicated. It is
sophisticated. It should require a degree of time, a certain amount of effort, and a certain amount
of motivation on behalf of the reader to find and discover insights.
[00:04:51.55] But my main three principles for visualisation-- it should be trustworthy, so the
audience should have a sense that what they are consuming in terms of information,
understanding, and knowledge is reliable. It is accurate. It is based on good judgement.
[00:05:09.98] The second principle is it should be accessible. And this does unpack the notion of
clarifying because accessible design is removing unnecessary obstacles to understanding. Now
again, that does not mean it's always about making things quick, simple, and easy. There is a
time and place for you to judge that that's important. But accessible design is about removing
unnecessary confusion, things that people don't understand how to read. What do those colours
represent? What do I do to click on that? How do I use this thing?
[00:05:43.81] The third principle in some respects is a bonus. It's about elegance. Can we make
our work as aesthetically appealing and as attractive and as seductive as possible? Now, one of
the interesting things about elegance as a concept in design is, you notice it when it's missing
more than when it's present. We take elegance for granted. So elegance is about the harmony of
the colour choices, about the layout, about the beauty of the typeface choices, things that will
cause someone to stop and to look at your work but will also help to sustain that engagement
throughout the process of them reading and understanding the work. So trustworthiness,
accessible, and elegance are the three things I think are most important.
[00:06:32.31] Oh, that's really, really important. And I think it elegantly summarises what we
should think of why doing good data visualisation. And I feel that-- when you're speaking, I feel
that there is a scientific perspective and a non-scientific one to make good data visualisation. So
from both these perspectives, do you think that when we want to present efficiently a particular
information we need to first understand how the human brain captures, understands, and
integrates efficiently a visual message?
[00:07:13.38] Absolutely. And I think the key word there is "human." And again, there are two
branches that we need to think about here because the human cognition process, for which there
is scientific evidence about things that work well and things that don't work well or at least
things that exploit the capability of the eye and the brain and the visual processing system. And
there are lots of studies through the ancestry of the field that tell us the things that we should do
and the things that we shouldn't do. 
[00:07:44.55] But the second branch that comes off the notion of "human" is people. People are
always the recipients of this. And people are inherently complex and flawed. We have blind
spots. We have prejudices. We have pressures. We have lack of confidence. We've already made
our mind up about something.
[00:08:04.86] So one of the absolute key fundamental bits of advice that any visualisation person
would always tell somebody else is, think about the audience. Think about the receiver. And
when we talk about that we're asking, what are the characteristics of these people? What do they
need to know? What do they currently not know? What are their motivations for what they will
do with this?
[00:08:28.83] Is it about direct decision making or actions? Or is it just an extra grain of
knowledge about something? Are they in a situation that's under pressure? It's noisy, it's busy.
They're travelling. They can't look at something in detail on the phones.
[00:08:43.20] So you're always trying to put yourself into the mindset of the recipient and think
about their capabilities and how they will encounter your work to understand some of these
decisions that you need to make. So it is absolutely a blend of the science but also the art. And by
art, I don't mean necessarily just creativity, but the art of judging these imperfections about
human responses and the situations that they go through. 
